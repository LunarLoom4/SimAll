"""hello.py — Smallest possible SimAll script.

Prints to the embedded Python Console panel and writes a value into the
shared ScriptContext bag.
"""
import simall

simall.log("Hello from the SimAll Python interpreter.")
ctx = simall.ctx()
ctx.set("examples.hello.ran", "true")
print("Context keys:", ctx.keys())
