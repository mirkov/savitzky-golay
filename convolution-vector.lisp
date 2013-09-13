;; Copyright Mirko Vukovic 2013.
;; Distributed under the Boost Software License, Version 1.0.
;; (See accompanying file ./LICENSE_1_0.txt or copy at
;; http://www.boost.org/LICENSE_1_0.txt)

;;;; convolution-vector.asd

(in-package :savitzky-golay)

(defun make-convolution-vector (n c-coeffs n-c n-r)
  "Store Savitzky-Golay coefficients C-COEFFS in a foreign vector of
length N, in a form suitable for convolution.  Also return (max n-r
n-l), where n-l is calculated as: n-c - n-r - 1

The Savitzky-Golay coefficients are of length N-C with N-R right
points (these are values provided by C-COEFFS).

It is assumed that c-coeffs is reversed, i.e., starting from point
n-r, n-r -1, n-r -2,...,0,-1,-2,...,- n-l

This vector can be FFT-ed and used in convolutions for data filtering

Coefficients n-r, n-r -1, n-r -2,..., are stored at the high-end of
the vector.  Coefficients  0,-1,-2,...,-n-l are stored at the low end of the vector"
  (let ((vec (make-foreign-array 'double-float
				 :dimensions n
				 :initial-element 0d0)))
    ;; store C-COEFFS in wrap-around order, starting at vector element
    ;; N - N-R.
    ;;
    ;; For example, if N-R=1, we store the first coefficient in the
    ;; last element of VEC, i.e., at N - 1
    (loop :for i :from (- n-r) :below (- n-c n-r)
       :do (setf (aref vec (mod (+ n i) n)) (pop c-coeffs)))
    (values vec (max n-r (- n-c n-r 1)))))

(defun convolution-vector (n p m k &optional (s 0))
  "Return a Savitzky-Golay filter for polynomial of order k, N-L left
points, N-R right points and derivative order S, in a vector of
length L-P.

The filter is stored in wrap-around order, ready to be FFT'd for
convolution purposes.

In wrap-around order, filter coefficient indices 0,...,N-L are stored
at the end of the vector in vector locations N - 1,...,N - N-L -
1.  Filter coefficients - N-L,...,-1 are stored at "
  (multiple-value-bind (c-coeffs n-c n-r)
      (filter-coeffs p m k s t)
    (make-convolution-vector n
			     (if (eql *method* :gram)
				 (mapcar (lambda (rational)
					    (float rational 1d0))
					  c-coeffs)
				 c-coeffs)
			     n-c n-r)))

(define-test convolution-vector-method-comparison
  "Compare the convolution result created by the two methods"
  (let ((*epsilon* 1e-5))
    (let ((m-vec (let ((*method* :matrix))
		   (copy-to (convolution-vector 10 0 2 2) 'array)))
	  (g-vec (let ((*method* :gram))
		   (copy-to (convolution-vector 10 0 2 2) 'array))))
      (assert-numerical-equal m-vec g-vec))
    (let ((m-vec (let ((*method* :matrix))
		   (copy-to (convolution-vector 10 1 2 2) 'array)))
	  (g-vec (let ((*method* :gram))
		   (copy-to (convolution-vector 10 1 2 2) 'array))))
      (assert-numerical-equal m-vec g-vec))
    (let ((m-vec (let ((*method* :matrix))
		   (copy-to (convolution-vector 10 1 2 2 1) 'array)))
	  (g-vec (let ((*method* :gram))
		   (copy-to (convolution-vector 10 1 2 2 1) 'array))))
      (assert-numerical-equal m-vec g-vec))))
