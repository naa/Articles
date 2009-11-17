;(define map-n map)
(define (zero? v)
	(= v 0))

(define (add a b)
  (map-n + a b))

(define (sum list)
  (fold-left add (car list) (cdr list)))

(define (sub a b)
  (map-n - a b))

(define (mul a b)
  (cond ((and (list? a) (number? b)) (map (lambda (x) (* x b)) a))
	((and (number? a) (list? b)) (mul b a))
	((and (number? a) (number? b)) (* a b))
	(else (error "Can not multiply"))))

(define (div a b)
  (map (lambda (x) (/ x b)) a))
  
(define (prod a b)
  (fold-left + 0 (map-n * a b)))

(define (simple-roots series rank)
  (cond ((eq? series 'A)
	 (do ((i 0 (+ i 1))
	      (base '()
		    (cons 
		     (do ((j 0 (+ j 1))
			  (v '()
			     (cons
			      (cond ((= j i) 1)
				    ((= j (+ i 1)) -1)
				    (else 0))
			      v)))
			 ((> j rank) (reverse v)))
		     base)))
	     ((= i rank) (reverse base))))
	((eq? series 'B)
	 (do ((i 0 (+ i 1))
	      (base '()
		    (cons 
		     (do ((j 0 (+ j 1))
			  (v '()
			     (cons
			      (cond ((= j i) 1)
				    ((= j (+ i 1)) -1)
				    (else 0))
			      v)))
			 ((= j rank) (reverse v)))
		     base)))
	     ((= i rank) (reverse base))))
	((eq? series 'C)
	 (do ((i 0 (+ i 1))
	      (base '()
		    (cons 
		     (do ((j 0 (+ j 1))
			  (v '()
			     (cons
			      (cond ((and (= i (- rank 1)) (= i j)) 2)
				    ((= j i) 1)
				    ((= j (+ i 1)) -1)
				    (else 0))
			      v)))
			 ((= j rank) (reverse v)))
		     base)))
	     ((= i rank) (reverse base))))
	((eq? series 'D)
	 (do ((i 0 (+ i 1))
	      (base '()
		    (cons 
		     (do ((j 0 (+ j 1))
			  (v '()
			     (cons
			      (cond ((and (= i (- rank 1)) (= (- i 1) j)) 1)
				    ((= j i) 1)
				    ((= j (+ i 1)) -1)
				    (else 0))
			      v)))
			 ((= j rank) (reverse v)))
		     base)))
	     ((= i rank) (reverse base))))))

(define (simple-co-roots series rank)
  (co-roots (simple-roots series rank)))

(define (co-root r)
  (div (mul 2 r) (prod r r)))

(define (co-roots roots)
  (map co-root roots))

(define (cartan-matrix simple-roots)
 (let ((cr (co-roots simple-roots)))
  (map (lambda (a)
	 (map (lambda (av) (prod a av)) cr))
	 simple-roots)))


(define (reflect weight root)
  (sub weight (mul (prod weight (co-root root)) root)))


(define (matrix-cofactor mat i j)
  (define (butnth n lst)
    (if (<= n 1)
	(cdr lst)
	(cons (car lst)
	      (butnth (+ -1 n) (cdr lst)))))
  (define (minor matrix i j)
    (map (lambda (x)
	   (butnth j x))
	 (butnth i mat)))
  (* (if (odd? (+ i j)) -1 1)
     (matrix-determinant (minor mat i j))))

(define (matrix-determinant mat)
  (let ((n (length mat)))
    (if (eqv? 1 n) (caar mat)
	(do ((j n (+ -1 j))
	     (ans 0 (+ ans (* (list-ref (car mat) (+ -1 j))
			      (matrix-cofactor mat 1 j)))))
	    ((<= j 0) ans)))))


(define (matrix-inverse mat)
  (let* ((det (matrix-determinant mat))
	 (rank (length mat)))
    (and (not (zero? det))
	 (do ((i rank (+ -1 i))
	      (inv '() (cons
			(do ((j rank (+ -1 j))
			     (row '()
				  (cons (/ (matrix-cofactor mat j i) det) row)))
			    ((<= j 0) row))
			inv)))
	     ((<= i 0) inv)))))

(define (fundamental-weights simple-roots)
  (map (lambda (x)
	 (sum (map-n (lambda (y z) (mul y z)) x simple-roots)))
       (matrix-inverse (cartan-matrix simple-roots))))

(define (rho simple-roots)
  (sum (fundamental-weights simple-roots)))

; orbit of Weyl group

(define (orbit weights roots)
  (let ((addon 
	 (filter (lambda (x) (not (element-of-set? x weights)))
		 (fold-right union-set '() 
			     (map (lambda (w)
				    (map (lambda (x) (reflect w x)) roots))
				  weights)))))
    (if (null? addon)
	weights
	(orbit (union-set weights addon) roots))))

(define (reflect-with-eps weight root)
  (if (or (null? (cdr weight))
	  (pair? (cdr weight)))
      (let ((v (reflect weight root)))
	(cons v (if (equal? weight v) 1 -1)))
      (let ((v (reflect (car weight) root)))
	(cons v (if (equal? (car weight) v) (cdr weight) (- (cdr weight)))))))

(define (orbit-with-eps weights roots)
  (let ((addon 
	 (filter (lambda (x) (not (element-of-set? x weights)))
		 (fold-right union-set '() 
			     (map (lambda (w)
				    (map (lambda (x) (reflect-with-eps w x)) roots))
				  weights)))))
    (if (null? addon)
	weights
	(orbit-with-eps (union-set weights addon) roots))))

(define (weight-orbit weight roots)
  (orbit-with-eps (list (cons weight 1)) roots))
  
(define (pos-roots simple-roots)
	(let ((r (rho simple-roots)))
		(filter (lambda (x) (>= (prod x r) 0)) 
			(orbit simple-roots simple-roots))))

(define (make-lie-algebra simple-roots)
  (lambda (tag)
    (cond
     ((equal? tag 'simple-roots) simple-roots)
     ((equal? tag 'rho)
      (let ((r (rho simple-roots)))
	r))
     ((equal? tag 'fundamental-weights)
      (let ((fw (fundamental-weights simple-roots))) 
	fw))
     ((equal? tag 'orbit)
      (lambda (weight) (orbit (list weight) simple-roots)))
     ((equal? tag 'orbit-with-eps)
      (lambda (weight) (orbit-with-eps (list (cons weight 1)) simple-roots))))))

(define (make-highest-weight-module highest-weight lie-algebra)
  (lambda (tag)
    (cond 
     ((equal? tag 'anomalous-points)
      (let ((ap
	     (map
	      (lambda (x) (cons (sub (car x) (lie-algebra 'rho)) (cdr x)))
	      ((lie-algebra 'orbit-with-eps)
	       (add highest-weight (lie-algebra 'rho))))))
	ap))
     ((equal? tag 'star)
      (let ((st
	     (map (lambda (x) (cons (sub (car x) (lie-algebra 'rho)) (cdr x)))
		  ((lie-algebra 'orbit-with-eps)
		   (lie-algebra 'rho)))))
	st)))))
      

;(define (anomalous-points 