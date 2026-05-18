# Fonts

Place the Inter font family here as:

```
fonts/inter/Inter-Regular.ttf
fonts/inter/Inter-Medium.ttf
fonts/inter/Inter-SemiBold.ttf
fonts/inter/Inter-Bold.ttf
```

Inter is released under the SIL Open Font License 1.1:
https://github.com/rsms/inter

`scripts/bootstrap_fonts.ps1` will download and unpack the official
release tarball into this directory on first build.

Segoe UI is shipped with Windows and used as the fallback. Inter +
Segoe UI together cover all UI text in the dark theme (spec §3.7).
