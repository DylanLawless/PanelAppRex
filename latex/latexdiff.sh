==== latex/latexdiff.sh ====
#!/bin/bash

shopt -s nullglob
rm -f ./qvApp2025lawless_marked_changes.*


# handle the section headings
latexdiff \
  --append-textcmd=section \
  --append-textcmd=subsection \
  --append-textcmd=subsubsection \
  --append-textcmd=paragraph \
  ./qvApp2025lawless_20250625_v0.tex \
  ./qvApp2025lawless.tex \
  > qvApp2025lawless_marked_changes.tex



