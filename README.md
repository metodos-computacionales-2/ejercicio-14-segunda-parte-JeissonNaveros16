Punto 1
Se obtiene una ecuación diferencial lineal de segundo orden.

Punto 2
La solución es de forma sinusoidal dado que el polinomio característico es: P(r)=mr^2+k, cuyas raices son complejas: r=+-sqrt(k/m)i. Luego la solución es de la forma: Acos(sqrt(k/m)t)+Bsin(sqrt(k/m)t), y como las condiciones iniciales son x(0)=1 y x'(0)=0 se obtiene: x(t)=cos(sqrt(k/m)t)

Punto 4
Como se puede ver en la gráfica "SolucionEDO.png" (a la izquierda) el método de euler no es el más indicado dado que acumula errores en la velocidad en cada iteración y la grafica de posición crece inesperadamente.

Punto 5
Ahora, en la misma figura "SolucionEDO.png" (a la izquierda) se ve que el método de Runge-Kutta se mantiene en constante oscilación sin crecer en magnitud, que es la solución esperada. Esto probablemente se debe a que antes de dar un avance para x, considera la contribución del avance en la velocidad (pues no es plenamente constante en el step), luego asocia pesos a sus 4 términos (k1, k2, k3 y k4) de avance.

Punto 6
En la figura "SolucionEDO.png" (a la derecha) se puede ver una elipse, esto se debe a que la posicion y la velocidad se comportan de manera sinusoidal: x(t)=cos(sqrt(k/m)t), x'(t)=-sqrt(k/m)sin(sqrt(k/m)t), de modo que ambas generan una elipse donde el semieje mayor está en el eje y (velocidad) debido al factor sqrt(k/m) (aprox 7,1). Por otro lado, el signo sólo indica que la elipse se recorre en sentido horario.

Punto 7
Cuando se incluye la fricción, se observa un comportamiento de la figura "SolucionEDOConFriccion.png" amortiguado. Esto se debe a que ahora la ecuación diferencial tiene polinomio característico: P(r)=mr^2+gr+k (g=gamma), cuyas raices tienen parte real: r=(-g+-sqrt(g^2-4km))/2m. Luego la solución general es de la forma: x(t)=exp(-g/2m)(Acos(sqrt(k/m-(g/2m)^2)t)+Bsin(sqrt(k/m-(g/2m)^2)t)). Finalmente se hace evidente que el factor exp(-g/2m) es el que porporciona el amortiguamiento.

Punto 8
Como se puede observar en la figura SolucionEDDistintosLambdas.png que dependiendo de lambda la función tiende a cambiar su frecuencia de oscilación (disminuyendo a mayor lambda).

Obs: Primero compilar y ejecutar el archipo .cpp (para generar los datos) y luego el .py (para generar las gráficas)