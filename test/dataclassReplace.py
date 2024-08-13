from dataclasses import dataclass, replace

@dataclass
class Foo:
  bar: str = "snafu"


def doit(f: Foo):
  print(f"doit: before {f=}")
  f = replace(f, bar = "foobar")
  print(f"doit: after {f=}")

foo = Foo()
print(f"before {foo=}")
doit(foo)
print(f"after {foo=}")
