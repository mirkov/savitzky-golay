;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; common.asd

(in-package :savitzky-golay)



(defparameter *method* :gram
"*METHOD* is used to set the default calculation method.  It can be
set to either :GRAM or :MATRIX.  Otherwise, the consequences are
undefined.")

(defgeneric weight-1 (method i p m k &optional s)
  (:documentation
"Calculate weight of the I-th data point for smoothing or
derivative calculation at the P-th point.  The weight is
calcualted for the S-th derivative of polynomial of order k over
2M+1 points.  

METHOD, :gram or :matrix, determines which algorithm is used to
calculate the weight.

All arguments are integers satisfying:
- |P| <= M
- |I| <= M
- M > 0
- K >= 0
- S >= 0.  The default is 0

The routine does not perform argument checks.
"))

(defun weight (i p m k &optional (s 0))
  "Calculate weight of the I-th data point for smoothing or
derivative calculation at the P-th point.  The weight is
calcualted for the S-th derivative of polynomial of order k over
2M+1 points.

The calculation method (Gram polynomial or matrix) is by the value
of *method*.

All arguments are integers satisfying:
- |P| <= M
- |I| <= M
- M > 0
- K >= 0
- S >= 0.  The default is 0

      
The routine does not perform argument checks.
"
  (weight-1 *method* i p m k s))

(defgeneric filter-coeffs-1 (method p m k &optional s reverse-p)
  (:documentation
"Returns a list of 2M+1 filter coefficients for smoothing or
derivative calculation at the P-th point.  The coefficients are
calculated for the S-th derivative of polynomial of order k over
2M+1 points.

METHOD, either +gram+ or +matrix+, determines which algorithm is
used to calculate the weight.

If /reverse-p/ is T, the filters are returned in reverse order

All arguments are integers satisfying:
- |P| <= M
- M > 0
- K >= 0
- S >= 0.  The default is 0

  
The routine does not perform argument checks.
"))

(defun filter-coeffs (p m k &optional s reverse-p)
"Returns a list of 2M+1 filter coefficients for smoothing or
derivative calculation at the P-th point.  The coefficients are
calculated for the S-th derivative of polynomial of order k over
2M+1 points.

The calculation method (Gram polynomial or matrix) is by the value
of *method*.

All arguments are integers satisfying:
- |P| <= M
- M > 0
- K >= 0
- S >= 0.  The default is 0

If REVERSE-P is T, the coefficients are returned in reverse order
      
The routine does not perform argument checks."

  (filter-coeffs-1 *method* p m k s reverse-p))
