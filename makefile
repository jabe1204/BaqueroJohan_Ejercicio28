ejercicio28.png : ejercicio28.dat graficador.py
	python graficador.py

ejercicio28.dat : ejercicio28.x
	./ejercicio28.x > ejercicio28.dat

ejercicio28.x : ejercicio28.cpp
	c++ ejercicio28.cpp -o ejercicio28.x
	
clean :
	rm ejercicio28.dat ejercicio28.x ejercicio28.png