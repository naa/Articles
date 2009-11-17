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