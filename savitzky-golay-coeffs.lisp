;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; savitzky-golay-coeffs.lisp

(in-package #:savitzky-golay)

(defun make-design-matrix (n-l n-r M)
  "Build design matrix, (NR-14.8.2)"
  (let* ((row-count (+ 1 n-l n-r))
	 (column-count (+ 1 M)) 
	 (A (make-foreign-array 'double-float
			       :dimensions (list row-count
						 column-count))))
    (loop
       :for i-row :from 0
       :for i :from (- n-l) :upto n-r
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
  "Basic test of coefficient calculations for n-l = n-r = m = 2

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

(defun c-n (n n-l n-r m &optional (l-d 0))
  "Return n-th coefficient of Savitzky-Golay filter for n-l, n-r and M
and differential order l-d

- n-l <= n <= n-r

This routine is used to test my implementation, and not in the
production code below

Unlike numerical recipes, I do not compute the matrix inverse.
Instead I solve for c-n rewriting NR-14.8.6 as

 (A^T . A ) c-n = A^T . e-n

I use LU decomposition to solve for c-n.  The result is 
"
  (let* ((A (make-design-matrix n-l n-r m))
	 (A^T-A (A^T-A A))
	 (n-absolute (+ n n-l))
	 (A^T-e-n (A^T-b A (make-unit-vector (+ 1 n-l n-r) n-absolute)))
	 (c-vec (multiple-value-bind (sign permutation)
		     (lu-decomposition A^T-A)
		   (declare (ignore sign))
		   (lu-solve A^T-A A^T-e-n permutation))))
    (aref c-vec l-d)))

(define-test c-coeffs
  "Test calculated coefficients against values of table on p.651

The assertions fail, as I don't understand how to control the tests.

I tried using the :test #'float-equal specification, but that did not
seem to affect the test"
  (let ((*significant-figures* 3))
    (assert-numerical-equal
     '(-0.086 0.343 0.486 0.343 -0.086)
     (c-coeffs 2 2 2))
    (assert-numerical-equal
     '(-0.143 0.171 0.343 0.371 0.257)
     (c-coeffs 3 1 2))
    (assert-numerical-equal
     '(0.086 -0.143 -0.086 0.257 0.886)
     (c-coeffs 4 0 2))
    (assert-numerical-equal
     '(-0.084 0.021 0.103 0.161 0.196 0.207 0.197 0.161 0.103 0.021 -0.084)
     (c-coeffs 5 5 2))
    (assert-numerical-equal
     '(0.035 -0.128 0.070 0.315 0.417 0.315 0.07 -0.128 0.035)
     (c-coeffs 4 4 4))
    (assert-numerical-equal
     '(0.042 -0.105 -0.023 0.140 0.280 0.333 0.280 0.14 -0.023 0.105 0.042)
     (c-coeffs 5 5 4))))

(defun c-coeffs (n-l n-r m &optional (l-d 0) (reverse t))
  "Return a list of Savitzky-Golay coefficients of polynomial order M,
N-L left points, N-R right points, and derivative order L-D.  Also
return number of coefficiens and n-r.

It returns a list of coefficients, starting with index
 - N-L,...,0,...,N-R

If REVERSE is T, the coefficients are reversed.  This is useful (and
required) for building convolution vectors"
  (assert (>= 0 n-l) ()
	  (error "n-l, ~a, must be zero or positive" n-l))
  (assert (>= 0 n-r) ()
	  (error "n-r, ~a, must be zero or positive" n-r))
  (assert (>= 0 m) ()
	  (error "m, ~a, must be zero or positive" m))
  (assert (>= 0 l-d) ()
	  (error "l-d, ~a, must be zero or positive" l-d))
  (assert (<= l-d m) ()
	  (error "l-d, ~a, must be equal to or smaller than m, ~a" l-d m))
  (assert (>= (+ n-l n-r) m) ()
	  (error "Sum of n-l, ~a, and n-r, ~a, must equal or exceed m, ~a"
		 n-l n-r m))
  (let* ((A (make-design-matrix n-l n-r m))
	 (A^T-A (A^T-A A))
	 (1+n-l+n-r (+ 1 n-l n-r))
	 c-coeffs)
    (multiple-value-bind (sign permutation)
	(lu-decomposition A^T-A)
      (declare (ignore sign))
      (dotimes (n 1+n-l+n-r c-coeffs)
	(let* ((A^T-e-n (A^T-b A (make-unit-vector 1+n-l+n-r n)))
	       (result 
		(lu-solve A^T-A A^T-e-n permutation)))
	  (push (aref result l-d) c-coeffs))))
    (values (if reverse
		c-coeffs
		(nreverse c-coeffs))
	    1+n-l+n-r n-r)))

(defun make-convolution-vector (n-p c-coeffs n-c n-r)
  "Store Savitzky-Golay coefficients C-COEFFS in a foreign vector of
length N-P, in a form suitable for convolution.

The Savitzky-Golay coefficients are of length N-C with N-R right
points (these are values provided by C-COEFFS).

It is assumed that c-coeffs is reversed, i.e., starting from point
n-r, n-r -1, n-r -2,...,0,-1,-2,...,- n-l

This vector can be FFT-ed and used in convolutions for data filtering

Coefficients n-r, n-r -1, n-r -2,..., are stored at the high-end of
the vector.  Coefficients  0,-1,-2,...,-n-l are stored at the low end of the vector"
  (let ((vec (make-foreign-array 'double-float
				 :dimensions n-p
				 :initial-element 0d0)))
    ;; store C-COEFFS in wrap-around order, starting at vector element
    ;; N-P - N-R.
    ;;
    ;; For example, if N-R=1, we store the first coefficient in the
    ;; last element of VEC, i.e., at N-P - 1
    (loop :for i :from (- n-r) :below (- n-c n-r)
       :do (setf (aref vec (mod (+ n-p i) n-p)) (pop c-coeffs)))
    vec))

(defun convolution-vector (n-p n-l n-r m &optional (l-d 0))
  "Return a Savitzky-Golay filter for polynomial of order M, N-L left
points, N-R right points and derivative order L-D, in a vector of
length L-P.

The filter is stored in wrap-around order, ready to be FFT'd for
convolution purposes.

In wrap-around order, filter coefficient indices 0,...,N-L are stored
at the end of the vector in vector locations N-P - 1,...,N-P - N-L -
1.  Filter coefficients - N-L,...,-1 are stored at "
  (multiple-value-bind (c-coeffs n-c n-r)
      (c-coeffs n-l n-r m l-d)
    (make-convolution-vector n-p c-coeffs n-c n-r)))
