 ; LispMeObjects
 ; http://c2.com/cgi/wiki?LispMeObjects
 ; written by Don Wells
 ; Create a new class with (class name super '(slot value)... '(method function)).
 ; Always use 'object as the super
 ; class at the very least.
 ; a function used as a method 
 ; will take at least one argument 
 ; self, the object that originally
 ; received the method.
 ; Invoke a function by sending the 
 ; name and arguments to an
 ; object. (e.g. (send anObject 'add 'sum 10))
 ; where add is the method and sum and 10
 ; are arguments)
 ; Get the value of a slot by sending
 ; the slot's name.
 ; (e.g. (send anObject 'sum))
 ; Set the value of a slot by sending
 ; the set method defined on object.
 ; (e.g. (send anObject 'set 'sum 20))
 ; Always evaluate (clearClasses) before
 ; doing anything.

 ; an object is (superName (slotname value)... (methodName closure)...)
 ; a class is (className . object)

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