(TeX-add-style-hook "documentation"
 (function
  (lambda ()
    (LaTeX-add-labels
     "sec:aux")
    (TeX-run-style-hooks
     "latex2e"
     "art10"
     "article"
     "interface"
     "param"
     "schedule"))))

