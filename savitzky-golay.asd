;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; savitzky-golay.asd

(asdf:defsystem #:savitzky-golay
  :serial t
  :description "Describe savitzky-golay here"
  :author "Your Name <your.name@example.com>"
  :license "Specify license here"
  :depends-on (#:lisp-unit
               #:antik
               #:gsll)
  :components ((:file "savitzky-golay-package-def")
               (:file "savitzky-golay-coeffs")))

