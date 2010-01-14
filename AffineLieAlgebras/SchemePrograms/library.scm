(define nil ())
(define true #t)
(define false #f)

(define (zero? v)
	(= v 0))

(define (map op list)
	(if (null? list) ()
		(cons (op (car list)) (map op (cdr list)))))

(define (map-n op . lists) 
    (define (map-n0 op  lists)
	(if (or (null? lists) (null? (car lists))) '()
	(cons (apply op (map car lists))
		(map-n0 op  (map cdr lists)))))
    (map-n0 op lists))

(define (accumulate op initial sequence)
  (if (null? sequence)
      initial
      (op (car sequence)
          (accumulate op initial (cdr sequence)))))

(define (filter predicate sequence)
  (cond ((null? sequence) ())
        ((predicate (car sequence))
         (cons (car sequence)
               (filter predicate (cdr sequence))))
        (else (filter predicate (cdr sequence)))))

(define (fold-left op initial sequence)
  (define (iter result rest)
    (if (null? rest)
        result
        (iter (op result (car rest))
              (cdr rest))))
  (iter initial sequence))

(define (fold-right op init s) (accumulate op init s))
(define (element-of-set? x set)
  (cond ((null? set) false)
        ((equal? x (car set)) true)
        (else (element-of-set? x (cdr set)))))

(define (adjoin-set x set)
  (if (element-of-set? x set)
      set
      (cons x set)))

(define (intersection-set set1 set2)
  (cond ((or (null? set1) (null? set2)) '())
        ((element-of-set? (car set1) set2)        
         (cons (car set1)
               (intersection-set (cdr set1) set2)))
        (else (intersection-set (cdr set1) set2))))

(define (union-set set1 set2)
	(if (null? set1) set2
		(union-set (cdr set1)
			   (adjoin-set (car set1) set2))
	)
)

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


(define *classes* '())


 (define (clearClasses)
  (set! *classes* 
    `((object #f 
      (set ,setSlot)
      (super
         ,(lambda (self)
            (getClass (car self))))))))


 (define (setSlot self aSlotName aValue)
  (let ((slot (assoc aSlotName (cdr self))))
    (cond
      ((not slot) 
        (set-cdr! self 
          (cons 
            (list aSlotName  aValue) 
            (cdr self))))
      (else
        (set-car! (cdr slot) aValue))))
  aValue)


 (define (getClass aClass)
  (let ((class (assoc aClass *classes*)))
      (cond
        ((not class) #f)
        (else (cdr class)))))


 (define (class aName aSuperName . aDefinition)
  (set! *classes* 
    (cons 
      (cons aName (cons aSuperName aDefinition))
      *classes*))
  aName)


 (define (new aSuperName . args)
  (cons aSuperName args))


 (define (send anObject aMessage . args)
  (sendWithSelf anObject anObject aMessage args))


 (define (sendWithSelf self anObject aMessage args)
  (let 
    ((superName (car anObject))
     (slot (assoc aMessage (cdr anObject))))
    (cond
      (slot (valueOfSlot self slot args))
      ((not superName) #f)
      (else 
        (let ((superClass (getClass superName)))
          (cond
            ((not superClass) #f)
            (else
              (sendWithSelf self superClass aMessage args))))))))


 (define (valueOfSlot self theSlot args)
  (let ((value (cadr theSlot)))
    (cond
      ((procedure? value)
        (apply value (cons self args)))
      (else value))))


