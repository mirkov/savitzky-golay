;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; gram-matrix-tests.lisp

(in-package #:savitzky-golay)
;; We compare coefficient values of the two methods



(define-test gram/matrix-coeff-comp
  "Test values calculated by the two methods.  This is the most
stringent test that I am doing.

I test for a second order filter coefficients:
- center point
- offest point
- derivative at the offset point"
  (let ((*epsilon* 1e-15))
    (assert-numerical-equal
     (multiple-value-bind (coeffs &rest ignore)
	 (filter-coeffs-1 :matrix 0 2 2)
       (declare (ignore ignore))
       coeffs)
     (mapcar (lambda (rational)
	       (float rational 1.d0)) (filter-coeffs-1 :gram 0 2 2)))
    (assert-numerical-equal
     (multiple-value-bind (coeffs &rest ignore)
	 (filter-coeffs-1 :matrix -1 2 2)
       (declare (ignore ignore))
       coeffs)
     (mapcar (lambda (rational)
	       (float rational 1.d0)) (filter-coeffs-1 :gram -1 2 2)))
    (assert-numerical-equal
     (multiple-value-bind (coeffs &rest ignore)
	 (filter-coeffs-1 :matrix -1 2 2 1)
       (declare (ignore ignore))
       coeffs)
     (mapcar (lambda (rational)
	       (float rational 1.d0)) (filter-coeffs-1 :gram -1 2 2 1)))))




