;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; savitzky-golay-coeffs.lisp

(in-package #:savitzky-golay)

(defun make-design-matrix (l r M)
  "Build design matrix, (NR-14.8.2)"
  (let* ((row-count (+ 1 l r))
	 (column-count (+ 1 M)) 
	 (A (make-foreign-array 'double-float
			       :dimensions (list row-count
						 column-count))))
    (loop
       :for i-row :from 0
       :for i :from (- l) :upto r
       :do (loop for j-column :from 0 :upto M
		for j :from 0
		:do (setf (aref A i-row j-column)
			  (float (expt i j) 1d0))))
    A))



(defun make-unit-vector (M n &optional (value 1d0))
  "Return a vector of length M of zeroes except for N-th element which
stores VALUE"
  (let ((vector (make-foreign-array 'double-float
				    :dimensions (list M)
				    :initial-element 0d0)))
    (setf (aref vector n) value)
    vector))

(defun A^T-A (A &optional (inplacep nil))
  "Return matrix product of A^T and A.

If INPLACEP is T, the original matrix A is used.  Otherwise, it is
left undisturbed"
  (if inplacep
      (matrix-product A A nil 1d0 1d0 :trans :notrans)
      (let ((A (copy-to A 'grid:foreign-array)))
	(matrix-product A A nil 1d0 1d0 :trans :notrans))))

(defun A^T-b (A b)
  "Return product of matrix transpose of A and vector b

If INPLACEP is T, the original matrix A is used.  Otherwise, it is
left undisturbed"
  (matrix-product A b nil 1d0 1d0 :trans))

(define-test c-n-2-2-2
  "Basic test of coefficient calculations for l = r = m = 2

I compare with first row of table on (NR-p.651)"
  (let ((*significant-figures* 3))
    (assert-sigfig-equal 0.486
			 (c-n 0 2 2 2))
    (assert-sigfig-equal 0.343
			 (c-n -1 2 2 2)))
  (let ((*significant-figures* 2))
    (assert-sigfig-equal -0.086
			 (c-n -2 2 2 2)))
  (assert-float-equal (c-n 1 2 2 2)
		      (c-n -1 2 2 2))
  (assert-float-equal (c-n 2 2 2 2)
		      (c-n -2 2 2 2)))

(define-test c-n-misc
  "Spot check of c-n coefficients with values in table NR-p.651"
  (let ((*significant-figures* 3))
    (assert-sigfig-equal -0.143 (c-n -3 3 1 2))
    (assert-sigfig-equal 0.161 (c-n -2 5 5 2))
    (assert-sigfig-equal -0.128 (c-n 3 4 4 4))
    (assert-sigfig-equal 0.280 (c-n 1 5 5 4)))
  (let ((*significant-figures* 2))
    (assert-sigfig-equal 0.086 (c-n -4 4 0 2))))

(defun c-n (i l r k &optional (s 0))
  "Return i-th coefficient of Savitzky-Golay filter for l, r and k
and differential order s

- l <= i <= r

This routine is used to test my implementation, and not in the
production code below

Unlike numerical recipes, I do not compute the matrix inverse.
Instead I solve for c-n rewriting NR-14.8.6 as

 (A^T . A ) c-n = A^T . e-n

I use LU decomposition to solve for c-n.  The result is 
"
  (let* ((A (make-design-matrix l r k))
	 (A^T-A (A^T-A A))
	 (n-absolute (+ i l))
	 (A^T-e-n (A^T-b A (make-unit-vector (+ 1 l r) n-absolute)))
	 (c-vec (multiple-value-bind (sign permutation)
		     (lu-decomposition A^T-A)
		   (declare (ignore sign))
		   (lu-solve A^T-A A^T-e-n permutation))))
    (aref c-vec s)))

(define-test matrix-coeffs
  "Test calculated coefficients against values of table on p.651 of
  Numerical Recipes

One of the assertions fails.  There ought to be a test for a number of digits.

I tried using the :test #'float-equal specification, but that did not
seem to affect the test"
  (let ((*epsilon* 1e-2))
    (assert-numerical-equal
     '(-0.086 0.343 0.486 0.343 -0.086)
     (matrix-coeffs 2 2 2))
    (assert-numerical-equal
     '(-0.143 0.171 0.343 0.371 0.257)
     (matrix-coeffs 3 1 2 0 nil))
    (assert-numerical-equal
     '(0.086 -0.143 -0.086 0.257 0.886)
     (matrix-coeffs 4 0 2 0 nil))
    (assert-numerical-equal
     '(-0.084 0.021 0.103 0.161 0.196 0.207 0.197 0.161 0.103 0.021 -0.084)
     (matrix-coeffs 5 5 2))
    (assert-numerical-equal
     '(0.035 -0.128 0.070 0.315 0.417 0.315 0.07 -0.128 0.035)
     (matrix-coeffs 4 4 4))
    (assert-numerical-equal
     '(0.042 -0.105 -0.023 0.140 0.280 0.333 0.280 0.14 -0.023 0.105 0.042)
     (matrix-coeffs 5 5 4))))

(defmethod weight-1 ((method (eql :matrix)) i p m k &optional (s 0))
  (c-n (+ i m) (+ m p) (- m p) k s))

(defun matrix-coeffs (l r k &optional (s 0) (reverse t))
  "Return a list of Savitzky-Golay coefficients of polynomial order k,
L left points, R right points, and derivative order S.  Also
return number of coefficiens and r.

It returns a list of coefficients, starting with index
 - L,...,0,...,R

If REVERSE is T, the coefficients are reversed.  This is useful (and
required) for building convolution vectors"
  (assert (>= l 0) ()
	  (error "l, ~a, must be zero or positive" l))
  (assert (>= r 0) ()
	  (error "r, ~a, must be zero or positive" r))
  (assert (>= k 0) ()
	  (error "k, ~a, must be zero or positive" k))
  (assert (>= s 0) ()
	  (error "s, ~a, must be zero or positive" s))
  (assert (<= s k) ()
	  (error "s, ~a, must be equal to or smaller than m, ~a" s k))
  (assert (>= (+ l r) k) ()
	  (error "Sum of l, ~a, and r, ~a, must equal or exceed m, ~a"
		 l r k))
  (let* ((A (make-design-matrix l r k))
	 (A^T-A (A^T-A A))
	 (1+l+r (+ 1 l r))
	 matrix-coeffs)
    (multiple-value-bind (sign permutation)
	(lu-decomposition A^T-A)
      (declare (ignore sign))
      (dotimes (n 1+l+r)
	(let* ((A^T-e-n (A^T-b A (make-unit-vector 1+l+r n)))
	       (result 
		(lu-solve A^T-A A^T-e-n permutation)))
	  (push (aref result s) matrix-coeffs))))
    (values (if reverse
		matrix-coeffs
		(nreverse matrix-coeffs))
	    1+l+r r)))

#+unnecessary(defun matrix-coeffs-1 (l r k &optional (s 0) (reverse-p t))
  "Identical to MATRIX-COEFFS except that it returns only the
coefficients, and not the size information"
  (multiple-value-bind (coeffs coeff-count r)
      (matrix-coeffs l r k s reverse-p)
    (declare (ignore coeff-count r))
    coeffs))

(defmethod filter-coeffs-1 ((method (eql :matrix)) p m k &optional (s 0) reverse-p)
  (matrix-coeffs (+ m p) (- m p) k s reverse-p))


