;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; gram-poly.lisp


;; The functions in this file compute the coefficients using the
;; orthogonal Gram polynomials method discussed by Gurry.  This method
;; works for an odd number of points, and returns the coefficients as
;; ratationals.

(in-package :savitzky-golay)

(defun gram-poly (i m k &optional (s 0))
  "Evaluate Gram polynomial (or it's s-th derivative) of order K over 2M+1
points at point I.

All arguments are integers.  M, K and S must be positive integers.  M
must be non-zero.  I is in the range -M to M.

According to this definition, the gram polynomials that I calculate
are not properly normalized.  See:
http://en.wikipedia.org/wiki/Discrete_Chebyshev_polynomials

However they do lead to the correct Savitzky-Golay filter coefficients
when tested against Gurry's article and the matrix method of
calculating them.  (see unit tests in the rest of the file)
"
  ;; Code modeled after Gurry, Figure I
  (if (> k 0)
      (- (* (/ (* 2 (- (* 2 k) 1))
	       (* k (+ (* 2 m) (- k) 1)))
	    (+ (* i (gram-poly i m (- k 1) s))
	       (* s (gram-poly i m (- k 1) (- s 1)))))
	 (* (/ (* (- k 1) (+ (* 2 m) k))
	       (* k (+ (* 2 m) (- k) 1)))
	    (gram-poly i m (- k 2) s)))
      (if (and (zerop k)
	       (zerop s))
	  1 0)))




(defun gram-poly-norm^2 (p0 p1 2m+1)
  "From the wikipedia page the norm is defined the square root of the
  inner product, which itself is defined as sum of the squares/
  2m+1 (where 2m+1 is the total number of points

http://en.wikipedia.org/wiki/Discrete_Chebyshev_polynomials
"
  (/ (reduce #'+ (mapcar #'* p0 p1))
     2m+1))

(define-test gram-poly-norm^2
  "Test for orthogonality and sum for Gram polynomials of orders 0, 1,
and 2 over 5 points (m=2)

The tests are as follows:
- test values for k=0 to be unity
- test sum of values for k # 0 to be zero
- test norms of k0, k1, k2
- test orthogonality of orders 0&1, 0&2 and 1&2

Tests of norms of k1 and k2 the definition on
http://en.wikipedia.org/wiki/Discrete_Chebyshev_polynomials
"
  (let* ((m 2)
	 (2m+1 (+ (* 2 m) 1))
	 (k0 (loop :for i :from -2 :to 2
		:collect (gram-poly i m 0 0)))
	 (k1 (loop :for i :from -2 :to 2
		:collect (gram-poly i m 1 0)))
	 (k2 (loop :for i :from -2 :to 2
		:collect (gram-poly i m 2 0))))
    (assert-true (every (lambda (value)
			  (= value 1))
			k0)
		 "k0 constant")
    (assert-equal 1 (gram-poly-norm^2 k0 k0 2m+1) "k0 norm")
    (assert-true (zerop (reduce #'+ k1)) "k1 sum")
    (assert-equal 1 (gram-poly-norm^2 k1 k1 2m+1) "k1 norm, expected to fail")
    (assert-true (zerop (reduce #'+ k2)) "k2 sum")
    (assert-equal 1 (gram-poly-norm^2 k2 k2 2m+1) "k2 norm, expected to fail")
    (assert-true (zerop (reduce #'+ (mapcar #'* k0 k1))) "k0 ortho k1")
    (assert-true (zerop (reduce #'+ (mapcar #'* k0 k2))) "k0 ortho k2")
    (assert-true (zerop (reduce #'+ (mapcar #'* k1 k2))) "k1 ortho k2")))


(define-test gram-poly-deriv
  (let* ((m 2)
	 (k0 (loop :for i :from -2 :to 2
		:collect (gram-poly i m 0 1)))
	 (k1 (loop :for i :from -2 :to 2
		:collect (gram-poly i m 1 1)))
	 (k2 (loop :for i :from -2 :to 2
		:collect (gram-poly i m 2 1))))
    (assert-true (every #'zerop k0))
    (assert-true (every (lambda (value)
			  (= value 1/2))
			k1))))


(defun gen-fact (a b)
  "Product of a, a-1, ..., a-b+1"
  (let ((gf 1))
    (loop :for j :from (+ a (- b) 1) :to a
       :do (setf gf (* gf j)))
    gf))

(defmethod weight-1 ((method (eql :gram)) i p m k &optional (s 0))
  ;; Code follows Gurry, Figure I
  (loop :for l :from 0 :to k
     :summing (* (+ (* 2 l) 1)
		 (/ (gen-fact (* 2 m) l)
		    (gen-fact (+ (* 2 m) l 1) (+ l 1)))
		 (gram-poly i m l 0)
		 (gram-poly p m l s))))

(defmethod filter-coeffs-1 ((method (eql :gram)) p m k &optional (s 0) reverse-p)
  (let ((coeffs (loop for i from (- m) to m
		   :collect (weight-1 method i p m k s))))
    (values (if reverse-p (nreverse coeffs) coeffs)
	    (1+ (* 2 m))
	    (- m p))))


(define-test smoothing-weights/center-point
  "Test output for quadratic smoothing (Table I)"
  ;; test for center point
  (assert-equal 17/35 (weight-1 :gram 0 0 2 2))
  (assert-equal 12/35 (weight-1 :gram 1 0 2 2))
  (assert-equal -3/35 (weight-1 :gram 2 0 2 2)))

(define-test smoothing-weights/leading-point
  "test for left point (P=-2)"
  (assert-equal 31/35 (weight-1 :gram -2 -2 2 2))
  (assert-equal 9/35 (weight-1 :gram -1 -2 2 2))
  (assert-equal -3/35 (weight-1 :gram 0 -2 2 2))
  (assert-equal -5/35 (weight-1 :gram 1 -2 2 2))
  (assert-equal 3/35 (weight-1 :gram 2 -2 2 2)))
(define-test smoothing-weights/symmetry
  (assert-equal (weight-1 :gram -2 -2 2 2) (weight-1 :gram 2 2 2 2))
  (assert-equal (weight-1 :gram -2 2 2 2) (weight-1 :gram 2 -2 2 2))
  (assert-equal (weight-1 :gram -2 0 2 2) (weight-1 :gram 2 0 2 2))
  (assert-equal (weight-1 :gram -1 0 2 2) (weight-1 :gram 1 0 2 2))
  (assert-equal (weight-1 :gram 0 -2 2 2) (weight-1 :gram 0 2 2 2)))

(define-test derivative-weights/center-point
  "test for center point"
  (assert-equal 0 (weight-1 :gram 0 0 2 2 1))
  (assert-equal -1/10 (weight-1 :gram -1 0 2 2 1))
  (assert-equal -2/10 (weight-1 :gram -2 0 2 2 1)))

(define-test derivative-weights/leading-point
  "test for center point"
  (assert-equal 40/70 (weight-1 :gram 0 -2 2 2 1))
  (assert-equal 13/70 (weight-1 :gram -1 -2 2 2 1))
  (assert-equal -54/70 (weight-1 :gram -2 -2 2 2 1)))

(define-test derivative-weights/symmetry
  (assert-equal (- (weight-1 :gram -2 0 2 2 1)) (weight-1 :gram 2 0 2 2 1))
  (assert-equal (- (weight-1 :gram 0 -2 2 2 1)) (weight-1 :gram 0 2 2 2 1))
  (assert-equal (- (weight-1 :gram -2 -2 2 2 1)) (weight-1 :gram 2 2 2 2 1))
  (assert-equal (- (weight-1 :gram -2 2 2 2 1)) (weight-1 :gram 2 -2 2 2 1)))

