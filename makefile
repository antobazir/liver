liver : liver.o initialize.o diff_adv_reac.o
	gcc -Ofast -D GEOM=2 -o liver -fopenmp -g liver.o initialize.o diff_adv_reac.o
	
diff_adv_reac.o: diff_adv_reac.c diff_adv_reac.h
	gcc -Ofast -o  diff_adv_reac.o -fopenmp -g -c diff_adv_reac.c -W -Wall -ansi -pedantic

initialize.o: initialize.c initialize.h
	gcc -Ofast -o  initialize.o -fopenmp -g -c initialize.c -W -Wall -ansi -pedantic	

liver.o : liver.c initialize.h
	gcc -Ofast -o  liver.o -fopenmp -g -c liver.c -W -Wall -ansi -pedantic

clean : 
	rm *.o
