#!/bin/bash

shopt -s nullglob
rm -f ./panelapprex2025lawless_marked_changes.*


# handle the section headings
latexdiff \
  --append-textcmd=section \
  --append-textcmd=subsection \
  --append-textcmd=subsubsection \
  --append-textcmd=paragraph \
  ./submission_version/panelapprex2025lawless_appnote_v0.tex \
  ./panelapprex2025lawless_appnote.tex \
  > panelapprex2025lawless_marked_changes.tex



