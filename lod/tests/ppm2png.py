#!/usr/bin/env python3
"""Convert P6 PPM frames in a sim output directory to PNG (stdlib only).

Usage: ppm2png.py <dir>
Writes <name>.png next to each .ppm. Uses zlib + manual PNG chunk encoding,
no third-party packages, so it runs under `python3 -I`.
"""
import os
import struct
import sys
import zlib


def read_ppm(path):
    with open(path, "rb") as f:
        data = f.read()
    if data[:2] != b"P6":
        raise ValueError("not a P6 PPM: " + path)
    i = 2
    vals = []
    while len(vals) < 3:
        while i < len(data) and data[i] in b" \t\r\n":
            i += 1
        if data[i:i + 1] == b"#":
            while i < len(data) and data[i] not in b"\r\n":
                i += 1
            continue
        j = i
        while j < len(data) and data[j] not in b" \t\r\n":
            j += 1
        vals.append(int(data[i:j]))
        i = j
    while i < len(data) and data[i] in b" \t\r\n":
        i += 1
    w, h, maxv = vals
    if maxv != 255:
        raise ValueError("only maxval 255 supported")
    px = data[i:i + w * h * 3]
    if len(px) != w * h * 3:
        raise ValueError("truncated pixel data: " + path)
    return w, h, px


def write_png(path, w, h, px):
    def chunk(tag, payload):
        out = struct.pack(">I", len(payload)) + tag + payload
        return out + struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF)

    raw = bytearray()
    stride = w * 3
    for y in range(h):
        raw.append(0)  # filter type 0 (none)
        raw += px[y * stride:(y + 1) * stride]
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)  # 8-bit RGB
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", ihdr))
        f.write(chunk(b"IDAT", zlib.compress(bytes(raw), 9)))
        f.write(chunk(b"IEND", b""))


def main(d):
    n = 0
    for name in sorted(os.listdir(d)):
        if not name.endswith(".ppm"):
            continue
        w, h, px = read_ppm(os.path.join(d, name))
        write_png(os.path.join(d, name[:-4] + ".png"), w, h, px)
        n += 1
    print("converted", n, "frames in", d)


if __name__ == "__main__":
    main(sys.argv[1])
