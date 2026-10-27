#set page(
  paper: "us-letter",
  margin: (x: 0.75in, y: 0.75in),
  numbering: "1 / 1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(90))
      #grid(
        columns: (1fr, auto),
        [*Erdős Problem \#1191 — Supplementary Verification Compendium*],
        [*DOI: 10.13140/RG.2.2.17420.27526*],
      )
      #v(-4pt)
      #line(length: 100%, stroke: 0.4pt + luma(160))
    ]
  }
)
#set text(font: "New Computer Modern", size: 10pt)
#set par(justify: true, leading: 0.58em)
#set heading(numbering: "1.")
#show raw: set text(font: "DejaVu Sans Mono", size: 7.8pt)
#show raw.where(block: true): block.with(
  fill: luma(248),
  stroke: 0.4pt + luma(200),
  inset: 8pt,
  radius: 3pt,
  width: 100%,
  breakable: true,
)

#align(center)[
  #text(size: 14.5pt, weight: "bold")[
    SUPPLEMENTARY VERIFICATION COMPENDIUM — ERDŐS PROBLEM \#1191
  ]
  #v(5pt)
  #text(size: 10.5pt, style: "italic")[
    Unconditional Diophantine Haar Closure, Exact Harmonic-Difference Conservation, and Multi-Octave Dyadic Budgets for Erdős Problem \#1191
  ]
  #v(5pt)
  *Author:* Lucas Parobek (`parobek@gmail.com`) #h(1em) | #h(1em) *Date:* October 27, 2026 #h(1em) | #h(1em) *Preprint DOI:* `10.13140/RG.2.2.17420.27526` \
  *Public Repository:* #link("https://github.com/doublegate/Erdos_1191")[`https://github.com/doublegate/Erdos_1191`]
]

#v(4pt)
#line(length: 100%, stroke: 0.6pt)
#v(4pt)

= Overview & Verification Instructions (`README_ancillary.txt`)

#raw(read("anc/README_ancillary.txt"), block: true, lang: "text")

= Lean 4 (`v4.29.0-rc8`) Zero-`sorry` Formal Verification Suite (`1191.lean`)

Run directly with `lake --wfail build` and `lake env lean Erdos1191/AxiomCheck.lean` (or `lean paper/anc/1191.lean`, 14 verified theorems, 0 errors, 0 warnings, 0 `sorry`s, 0 custom axioms):

#raw(read("anc/1191.lean"), block: true, lang: "lean")

= Standalone Python 3 Verification Suite (`verify_erdos1191_diophantine_haar.py`)

Run directly with `python3 paper/anc/verify_erdos1191_diophantine_haar.py` (Python 3 standard library only, 0 external dependencies):

#raw(read("anc/verify_erdos1191_diophantine_haar.py"), block: true, lang: "python")

= Machine-Readable Verification Certificate (`diophantine_haar_certificates.json`)

#raw(read("anc/diophantine_haar_certificates.json"), block: true, lang: "json")
