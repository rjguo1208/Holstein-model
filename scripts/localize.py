"""Build English HTML from shared page sources and a reviewed translation catalog.

Source text is the translation key: changing Chinese copy requires an updated
English entry. Markup, numerical tables, equations and data links stay shared.
Run with --extract to list the translatable strings grouped by source page.
"""
import argparse
from html import escape, unescape
from html.parser import HTMLParser
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
HAN = re.compile(r"[\u3400-\u9fff]")
ATTRIBUTE = re.compile(r'''([\w:-]+)(\s*=\s*)(["'])(.*?)\3''', re.DOTALL)
MATH = re.compile(r"\\\[([\s\S]*?)\\\]|\\\(([\s\S]*?)\\\)")
PUNCTUATION = str.maketrans({
    "。": ".", "，": ", ", "、": ", ", "；": "; ", "：": ": ",
    "（": " (", "）": ")", "“": '"', "”": '"', "‘": "'", "’": "'",
    "！": "!", "？": "?", "／": "/",
})


def key(text):
    return " ".join(text.split())


def math_signature(source):
    # Only explanatory text inside equations is localized; all mathematics
    # must remain identical, including signs, constants and delimiters.
    return sorted(re.sub(r"\s+", "", re.sub(r"\\text\{[^{}]*\}", r"\\text{}", display or inline))
                  for display, inline in MATH.findall(source))


class Localizer(HTMLParser):
    def __init__(self, catalog=None):
        super().__init__(convert_charrefs=False)
        self.catalog = catalog
        self.parts = []
        self.keys = []

    def feed(self, data):
        self.source_lines = data.splitlines(keepends=True)
        super().feed(data)

    def translate(self, text):
        if not HAN.search(text):
            return text.translate(PUNCTUATION)
        normalized = key(text)
        self.keys.append(normalized)
        if self.catalog is None:
            return text
        translated = self.catalog.get(normalized)
        if not isinstance(translated, str) or not translated.strip():
            raise ValueError("Missing English translation: " + normalized)
        prose = MATH.sub("", translated)
        if HAN.search(translated) or "<" in prose or ">" in prose:
            raise ValueError("English translation contains Chinese text or HTML: " + normalized)
        prefix = text[:len(text) - len(text.lstrip())]
        suffix = text[len(text.rstrip()):]
        return prefix + translated.strip() + suffix

    def start(self, tag):
        if tag in {"a", "strong", "em", "code", "b", "i"} and self.parts:
            last = self.parts[-1]
            if last and not last[-1].isspace() and not last.endswith((">", "(")):
                self.parts.append(" ")
        def attribute(match):
            name, equals, quote, raw = match.groups()
            value = unescape(raw)
            if tag == "html" and name == "lang":
                value = "en"
            elif name in {"alt", "title", "aria-label", "content"}:
                value = self.translate(value)
            elif name in {"href", "src"} and value.startswith(("assets/", "figures/", "results/", "data/")):
                value = "../" + value
            else:
                return match.group(0)
            return name + equals + quote + escape(value, quote=True) + quote
        self.parts.append(ATTRIBUTE.sub(attribute, self.get_starttag_text()))

    def handle_starttag(self, tag, attrs):
        self.start(tag)

    def handle_startendtag(self, tag, attrs):
        self.start(tag)

    def handle_endtag(self, tag):
        self.parts.append("</" + tag + ">")

    def handle_data(self, data):
        translated = self.translate(data)
        if (self.parts and re.search(r"</(?:a|strong|em|code|span|b|i)>$", self.parts[-1])
                and translated and not translated[0].isspace()
                and translated[0] not in ".,:;!?)]}"):
            translated = " " + translated
        self.parts.append(translated)

    def handle_decl(self, decl):
        self.parts.append("<!" + decl + ">")

    def handle_comment(self, data):
        self.parts.append("<!--" + data + "-->")

    def handle_entityref(self, name):
        self.reference("&" + name)

    def handle_charref(self, name):
        self.reference("&#" + name)

    def reference(self, token):
        # HTMLParser also reports unterminated references such as LaTeX's
        # alignment token '&g^2'. Preserve the original spelling exactly.
        line, column = self.getpos()
        suffix = ";" if self.source_lines[line - 1][column:].startswith(token + ";") else ""
        self.parts.append(token + suffix)


def build(catalog=None):
    pages, strings = {}, {}
    for source in sorted((ROOT / "src").glob("*.html")):
        original = source.read_text(encoding="utf-8")
        parser = Localizer(catalog)
        try:
            parser.feed(original)
            parser.close()
        except ValueError as error:
            raise ValueError(source.name + ": " + str(error)) from error
        translated = "".join(parser.parts)
        if catalog is not None and math_signature(original) != math_signature(translated):
            raise ValueError(source.name + ": translation changed a mathematical expression")
        pages[source.name] = translated
        strings[source.name] = list(dict.fromkeys(parser.keys))
    if catalog is not None:
        used = {item for group in strings.values() for item in group}
        unused = set(catalog) - used
        if unused:
            raise ValueError("Remove or update obsolete translations: " + "; ".join(sorted(unused)))
    return pages, strings


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--extract", action="store_true")
    args = cli.parse_args()
    if args.extract:
        _, result = build()
    else:
        catalog = json.loads((ROOT / "src/locales/en.json").read_text(encoding="utf-8"))
        result, _ = build(catalog)
    print(json.dumps(result, ensure_ascii=False, indent=2))
