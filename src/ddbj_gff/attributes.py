from __future__ import annotations

from urllib.parse import unquote

_RESERVED = {
    "%": "%25",  # must come logically first when encoding
    ";": "%3B",
    "=": "%3D",
    "&": "%26",
    ",": "%2C",
    "\t": "%09",
    "\n": "%0A",
    "\r": "%0D",
}


def encode_value(text: str) -> str:
    out: list[str] = []
    for ch in text:
        if ch in _RESERVED:
            out.append(_RESERVED[ch])
        elif ord(ch) < 0x20:
            out.append("%%%02X" % ord(ch))
        else:
            out.append(ch)
    return "".join(out)


def decode_value(text: str) -> str:
    return unquote(text)


def parse_attributes(col9: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = {}
    if col9 in ("", "."):
        return attrs
    for pair in col9.split(";"):
        if pair == "":
            continue
        if "=" not in pair:
            key = decode_value(pair)
            attrs.setdefault(key, [])
            continue
        raw_key, raw_val = pair.split("=", 1)
        key = decode_value(raw_key)
        values = [decode_value(v) for v in raw_val.split(",")]
        if key in attrs:
            attrs[key].extend(values)
        else:
            attrs[key] = values
    return attrs


def serialize_attributes(attrs: dict[str, list[str]]) -> str:
    parts: list[str] = []
    for key, values in attrs.items():
        ek = encode_value(key)
        ev = ",".join(encode_value(v) for v in values)
        parts.append(f"{ek}={ev}")
    return ";".join(parts)
