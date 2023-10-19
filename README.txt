==============
auto r8.py
==============

MEJORAS
(1) Se corrigió el problema de los indices al NO procesar los archivos con Elongaciones fuera del rango
Se definio S como Arreglo en lugar de Lista. 

Esto en las primeras versiones NO era un problema porque TODOS los archivos se procesaban.

(2) Se especificó los dPIs en el tamaño de las figuras. 
dpi=150 que es lo usado para 'ebooks'

(3) Se hicieron pequeñas modificaciones en los prints y se añadio la cantidad de archivos a procesar.

==============
auto r9.py
==============

MEJORAS
(1) Se corrigió el problema de que si NO encuentra la radiofuente en el archivo de "radio sources.txt"
    entonces se identifique con su Ascención Recta & Declinación, tanto en la gráfica como en el archivo WIPSS
(2) Cambios en la figura de la Potencia
    a. No restringuí a 7x5 pulgadas
	b. El máximo en lo reduje de 10 a 6, para que no quedará tan alta
	c. Movi algo las coordenadas de los datos a desplegar
	d. Cambie ligeramente la distribución de los espacios en el Título de la gráfica
	e. Fije las leyendas en la esquina superior izquierda, para evitar los empalmes con los datos
==============
auto r10.py
==============

MEJORAS
(1) Se agregó la validación de que si todas las radiofuentes se salen de rango NO haga el MAPA de radiofuentes
(2) Se añadió un RESUMEN al final sobre la cantidad de radiofuentes en rango, fuera de rango, con Nans y NO identificadas
(3) Además de desplegar a pantalla ahora de guarda en un archivo LOG 
    En ./wipss/logs/log-<hora local>.txt

==============
auto r11.py
==============

MEJORAS
(1) Para el AJUSTE de Minimos Cuadrados 
    a. Se definió incluir a la Razón Axial como parámetro (NOTA: bajar nuevamente FUNCTIONS)
	b. Las condiciones iniciales son:
		Velocidad = 285,000    Alpha = 3.5     AR = 1
	c. Los límites son:
		Velocidad = [100,000;2,000,000]
		Alpha = [3.3;3.8]
		AR = [0.7;1.3]
	d. Dejar fija el Ancho de la Radiofuente a 0.2
		width = 0.2
(2) Cambios en el formato WIPSS
	a. Se detectó faltante las últimas 6 columnas que son: Method Vel. V-err g-value g-err Method
	b. En Method se pondrá SS como abreviación de Single Station
	c. Cambio en el Ancho de la Radio Fuente de 0.25 a 0.2
	d. Se definió para las Radiofuentes NO identificadas seguir el estandard de la IAU
        Posición 0 -> J	(por tomar como base el 2000)
		Posiciones 1-4 -> Los primeros 4 números de la Ascensión Recta
		Posiciones 5-9 -> Los primeros 4 números de la Declinación, incluyendo signo
	e. Al NOMBRE DEL ARCHIVO se le añadió la HORA LOCAL y se eliminó Limpiar Carpetas
(3) Modificación al nombre otorgado al Mapa de las Radiofuentes
			test-<fecha>  ->  maps-<fecha>-<hora local>    La FECHA es cuando se observaron las radiofuentes 
(4) Se eliminó salvar en archivos Postscript
(5) Se corrigió instrucción mal identada que provocaba la selección equivocada del tamaño del 
	círculo que representa la radiofuente.
    (Error introducido al agregar el if de NO graficar el MAPA si todas estaba fuera del rango válido para la elongación).

==============
auto-mexart.py.py
==============
(1) Identico a auto r11.py con excepción de tres aspectos:
	a. El(los) archivo(s) a analizar se deben de pasar en la línea de comandos al ejecutar el programa
		Por ejemplo:
		#python auto-mexart.py 20220302_005408-033355.dat 20220302_013741+330935.dat 20220302_000623-000047.dat
		
		O dentro del Spyder
		#run auto-mexart.py 20220302_005408-033355.dat 20220302_013741+330935.dat 20220302_000623-000047.dat

(2) Al nombre del archivo Wipss, que es WIPSS-FORMAT-FILE no se le añade la hora local, es decir que se preserva 
	en la carpeta ./wipss/list un SOLO archivo a menos que el OPERADOR lo copie o salve con otro nombre
	a. En cada corrida se añadirán más entradas a dicho archivo.
(3) El MAPA de radiofuentes sólo contendrá las ingresadas en DICHA CORRIDA.


==============
NOTA
	Los achivos DESTACADOS son las últimas versiones
