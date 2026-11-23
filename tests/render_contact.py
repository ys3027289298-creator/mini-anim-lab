#!/usr/bin/env python3
"""Render the four pose.ppm files into one PNG contact sheet (stdlib only).

PPM P6 -> PNG via zlib, used purely for visual acceptance of the new
pose export images. Run with: python3 -I tests/render_contact.py <ppm...>
"""
import struct
import sys
import zlib
from pathlib import Path


def read_ppm(path):
    data = Path(path).read_bytes()
    assert data[:2] == b"P6"
    # tokenize the header (max 3 whitespace-separated tokens, rest is binary)
    idx = 2
    tokens = []
    while len(tokens) < 3:
        while data[idx] in b" \t\n\r":
            idx += 1
        start = idx
        while data[idx] not in b" \t\n\r":
            idx += 1
        tokens.append(int(data[start:idx]))
    idx += 1  # single whitespace after maxval
    w, h, maxval = tokens
    assert maxval == 255
    return w, h, data[idx:idx + w * h * 3]


def png_chunk(tag, payload):
    return (struct.pack(">I", len(payload)) + tag + payload +
            struct.pack(">I", zlib.crc32(tag + payload) & 0xFFFFFFFF))


def write_png(path, w, h, rgb):
    raw = bytearray()
    stride = w * 3
    for y in range(h):
        raw.append(0)  # filter type 0
        raw.extend(rgb[y * stride:(y + 1) * stride])
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)  # 8-bit RGB
    out = (b"\x89PNG\r\n\x1a\n" +
           png_chunk(b"IHDR", ihdr) +
           png_chunk(b"IDAT", zlib.compress(bytes(raw), 9)) +
           png_chunk(b"IEND", b""))
    Path(path).write_bytes(out)


def main():
    ppms = [read_ppm(p) for p in sys.argv[1:]]
    w, h = ppms[0][0], ppms[0][1]
    pad = 8
    sheet_w = w * 2 + pad * 3
    sheet_h = h * 2 + pad * 3
    sheet = bytearray(b"\xff" * sheet_w * sheet_h * 3)
    for i, (cw, ch, rgb) in enumerate(ppms):
        col, row = i % 2, i // 2
        ox, oy = pad + col * (w + pad), pad + row * (h + pad)
        for y in range(ch):
            dst = ((oy + y) * sheet_w + ox) * 3
            src = y * cw * 3
            sheet[dst:dst + cw * 3] = rgb[src:src + cw * 3]
    out = Path(__file__).resolve().parent.parent / "artifacts" / "contact_sheet.png"
    write_png(out, sheet_w, sheet_h, sheet)
    print(out)


if __name__ == "__main__":
    main()
