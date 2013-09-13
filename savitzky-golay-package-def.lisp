;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; package.lisp



(defpackage #:savitzky-golay
  (:nicknames :sg :savgol)
  (:use #:cl #:grid #:gsll #:lisp-unit)
  (:shadow :cl :aref)
  (:shadow :lisp-unit :norm)
  (:shadow :gsll :row :column)
  (:export :c-coeffs
	   :make-convolution-vector
	   :convolution-vector)
  (:documentation "Package for calculating Savitzky-Golay coefficients
and returning them as a vector suitable for convolution

The package symbol naming is such that it is best not to import the
package, but to use the package nickname SG.

Thus, (sg:convolution-vector 10 3 1 2)"))


(antik:make-user-package :savitzky-golay)
