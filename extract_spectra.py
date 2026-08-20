#!/usr/bin/env python3
# extract_spectra.py
# Extract raw UV-Vis spectra from an Agilent/HP ChemStation .SD or .KD file.
# Written to run on Python 3.2.
#
# Usage:
#     python extract_spectra.py <filename.SD|filename.KD>
#
# Output: writes "<same base name>.csv" next to the input file.
# Headers: .SD -> spectrum name (blank -> "unnamed"); .KD -> elapsed seconds (plain number).

import csv
import os
import struct
import re
import sys


def read_double_field(data, name, start, end):
    pat = name.encode("utf-16-le")
    idx = data.find(pat, start, end)
    if idx == -1:
        return None
    after = idx + len(pat)
    return struct.unpack("<d", data[after + 6:after + 14])[0]


def read_string_field(data, name, start, end):
    pat = name.encode("utf-16-le")
    idx = data.find(pat, start, end)
    if idx == -1:
        return None
    after = idx + len(pat)
    length = struct.unpack("<H", data[after + 6:after + 8])[0]
    return data[after + 8:after + 8 + length * 2].decode("utf-16-le", "replace")


def read_sliced_axis(data, label, start, end):
    pat = label.encode("utf-16-le")
    idx = data.find(pat, start, end)
    if idx == -1:
        return None
    after = idx + len(pat)
    count = struct.unpack("<I", data[after + 4:after + 8])[0]
    rest = data[after + 8:after + 8 + 2 + 32]
    vals = struct.unpack("<4d", rest[2:34])
    axis_start = vals[2]
    axis_step = vals[3]
    return {"count": count, "start": axis_start, "step": axis_step}


def read_double_array(data, label, start, end):
    pat = label.encode("utf-16-le")
    idx = data.find(pat, start, end)
    if idx == -1:
        return None
    after = idx + len(pat)
    count = struct.unpack("<I", data[after + 4:after + 8])[0]
    arr_start = after + 9
    return struct.unpack("<" + str(count) + "d", data[arr_start:arr_start + 8 * count])


def extract(path):
    f = open(path, "rb")
    try:
        data = f.read()
    finally:
        f.close()

    wl_axis = read_sliced_axis(data, "Wavelength (nm)", 0, len(data))
    if wl_axis is None:
        raise ValueError("Could not find a 'Wavelength (nm)' axis - "
                         "is this really a ChemStation .SD/.KD UV-Vis file?")
    wavelengths = []
    for i in range(wl_axis["count"]):
        wavelengths.append(wl_axis["start"] + wl_axis["step"] * i)

    starts = []
    for m in re.finditer(re.escape(b"CHPUVObject"), data):
        starts.append(m.start())
    if not starts:
        raise ValueError("No 'CHPUVObject' spectrum records found.")
    starts.append(len(data))

    scans = []
    for i in range(len(starts) - 1):
        s = starts[i]
        e = starts[i + 1]
        absorbance = read_double_array(data, "Absorbance (AU)", s, e)
        if absorbance is None:
            continue
        sample_name = read_string_field(data, "SampleName", s, e)
        rel_time = read_double_field(data, "RelTime", s, e)
        scans.append((sample_name, rel_time, absorbance))

    return wavelengths, scans


def build_headers(scans, ext):
    raw = []
    if ext == ".KD":
        for item in scans:
            rel_time = item[1]
            if rel_time is not None:
                raw.append("%g" % rel_time)
            else:
                raw.append("")
    else:
        for item in scans:
            name = item[0]
            if name:
                raw.append(name)
            else:
                raw.append("unnamed")

    seen = {}
    headers = []
    for label in raw:
        seen[label] = seen.get(label, 0) + 1
        if seen[label] == 1:
            headers.append(label)
        else:
            headers.append("%s_%d" % (label, seen[label]))
    return headers


def write_csv(out_path, wavelengths, scans, ext):
    headers = build_headers(scans, ext)
    f = open(out_path, "w", newline="")
    try:
        w = csv.writer(f)
        w.writerow(["Wavelength_nm"] + headers)
        for row_i in range(len(wavelengths)):
            row = ["%.0f" % wavelengths[row_i]]
            for item in scans:
                absorbance = item[2]
                row.append("%.5f" % absorbance[row_i])
            w.writerow(row)
    finally:
        f.close()


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python extract_spectra.py <filename.SD|filename.KD>\n")
        sys.exit(1)

    in_path = sys.argv[1]
    if not os.path.isfile(in_path):
        sys.stderr.write("ERROR: file not found: %s\n" % in_path)
        sys.exit(1)

    ext = os.path.splitext(in_path)[1].upper()
    if ext not in (".SD", ".KD"):
        sys.stderr.write("ERROR: unsupported extension '%s'. Expected .SD or .KD.\n" % ext)
        sys.exit(1)

    out_path = os.path.splitext(in_path)[0] + ".csv"

    print("Reading %s ..." % os.path.basename(in_path))
    wavelengths, scans = extract(in_path)
    print("Found %d spectra, %d wavelength points each (%.0f-%.0f nm)." % (
        len(scans), len(wavelengths), wavelengths[0], wavelengths[-1]))

    write_csv(out_path, wavelengths, scans, ext)
    print("Wrote %s" % out_path)


if __name__ == "__main__":
    main()