(TeX-add-style-hook "countreg"
 (lambda ()
    (LaTeX-add-bibliographies)
    (LaTeX-add-labels
     "sec:intro"
     "sec:software"
     "eq:family"
     "eq:mean"
     "eq:Poisson"
     "eq:negbin"
     "eq:zeroinfl"
     "eq:zeroinfl-mean"
     "eq:hurdle"
     "eq:hurdle-mean"
     "sec:illustrations"
     "fig:ofp"
     "fig:bad-good"
     "fig:ofp2"
     "sec:summary"
     "app:zeroinfl"
     "app:hurdle"
     "app:methods")
    (TeX-add-symbols
     '("fct" 1)
     '("class" 1))
    (TeX-run-style-hooks
     "thumbpdf"
     "latex2e"
     "Z10"
     "Z")))

