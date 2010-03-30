(class 'lie-algebra 'object
       `(rank ,(lambda (self)
       	      (error "Abstract class!")))
       `(dimension ,(lambda (self)
       	      (error "Abstract class!")))
       `(cartan-subalgebra ,(lambda (self) 
       	      (error "Abstract class!")))
)
(class 'commutative-algebra 'lie-algebra
       '(dim 1)
       `(rank ,(lambda (self) (send self 'dim)))
       `(dimension ,(lambda (self) (send self 'dim)))
       `(cartan-subalgebra ,(lambda (self) self)))

;; The constructor
(define (make-commutative-algebra dim)
  (new 'commutative-algebra `(dim ,dim)))
(class 'semisimple-lie-algebra 'lie-algebra
       '(simple-roots ())
       `(rank ,(lambda (self)
           (length (send self 'simple-roots))))
              
       `(cartan-subalgebra ,(lambda (self)
           (make-commutative-algebra (send self 'rank))))
       `(simple-co-roots ,(lambda (self)
       	      (map co-root (send self 'simple-roots))))
       `(cartan-matrix
         ,(lambda (self)
            (map (lambda (a)
       	    (map (lambda (av) (prod a av)) (send self 'simple-co-roots)))
       	 (send self 'simple-roots))))
       `(general-orbit , (lambda (self weights reflect)
           (let ((addon 
       	   (filter (lambda (x) (not (element-of-set? x weights)))
       		   (fold-right union-set '() 
       			       (map (lambda (w)
       				      (map
       				       (lambda (x) (reflect w x))
       				       (send self 'simple-roots)))
       				    weights)))))
             (if (null? addon)
       	  weights
       	  (send self 'general-orbit (union-set weights addon) reflect )))))

       `(orbit ,
         (lambda (self weights)
           (send self 'general-orbit weights reflect)))

       `(orbit-with-eps ,
         (lambda (self weights)
           (send self 'general-orbit weights reflect-with-eps)))
       `(fundamental-weights
         ,(lambda (self)
            (map (lambda (x)
       	    (sum (map-n (lambda (y z) (mul y z)) x (send self 'simple-roots))))
       	  (matrix-inverse (send self 'cartan-matrix)))))
       `(rho ,(lambda (self)
       	 (sum (send self 'fundamental-weights))))

       `(roots ,(lambda (self)
       	   (send self 'orbit (send self 'simple-roots)))))
;; The constructor
(define (make-semisimple-lie-algebra simple-roots)
  (new 'semisimple-lie-algebra `(simple-roots ,simple-roots)))

(define (co-root r)
  (div (mul 2 r) (prod r r)))

(define (reflect weight root)
  (sub weight (mul (prod weight (co-root root)) root)))

(define (reflect-with-eps weight root)
  (if (or (null? (cdr weight))
	  (pair? (cdr weight)))
      (let ((v (reflect weight root)))
	(cons v (if (equal? weight v) 1 -1)))
      (let ((v (reflect (car weight) root)))
	(cons v (if (equal? (car weight) v) (cdr weight) (- (cdr weight)))))))

(define (eps weight)
  (if (or (null? (cdr weight))
	  (pair? (cdr weight)))
      1
      (cdr weight)))

(define (no-eps weight)
  (if (or (null? (cdr weight))
	  (pair? (cdr weight)))
      weight
      (car weight)))
(define (make-simple-lie-algebra series rank)
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
  (make-semisimple-lie-algebra (simple-roots series rank)))

(class 'representation 'object
       `(dim ())
       `(multiplicity , (lambda (self weight)
			  (error "Class is abstract!")))
       '(lie-algebra ())
       )

(class 'highest-weight-representation 'representation
       '(highest-weight ())
              `(anomalous-weights , 
       	 (lambda (self)
       	   (let* ((algebra (send self 'lie-algebra))
       		  (rho (send algebra 'rho)))
       	     (map (lambda (x)
       		    (sub x rho))
       		  (send algebra 'orbit (add (send self 'highest-weight) rho))))))
       `(multiplicity , (lambda (self weight)
       		   (let ((star ))
       		     (-
       		      (sum (map (lambda (wg) (* (eps wg) (delta weight wg)))))
       		      (sum (map (lambda (wg)
       				    (* (eps wg)
       				       (send self 'multiplicity (sub weight wg))))
       				  star))))))
         
       	    )
