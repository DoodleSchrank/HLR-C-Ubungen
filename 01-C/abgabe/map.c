#include <stdio.h>

// Definieren Sie ein enum cardd
enum cardd { N=1, E=2, S=4, W=8 };

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
static char *map[ 3 ][ 3 ] =	{{"0", "0", "0" },
				{ "0", "0", "0" },
				{ "0", "0", "0" }};

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, enum cardd dir)
{
	// Alle Werte > 12 sind illegal, genauso wie die Werte 5,6,7 und 11.
	// Außerdem dürfen x und y nicht 2 überschreiten, da die Matrix nur von 0 bis 2 reicht.
	if(dir <= 12 && dir != 5 && dir != 6 && dir != 7 && dir != 11 && x < 3 && y < 3)
	{
		
		char *c = "0";
		switch(dir)
		{
			case 1: c = "N"; break;
			case 2: c = "E"; break;
			case 3: c = "NE"; break;
			case 4: c = "S"; break;
			case 8: c = " W"; break;
			case 9: c = "NW"; break;
			case 10: c = "SE"; break;
			case 12: c = "SW"; break;
			default: c = "0"; break;
		}
		map [ x ][ y ] = c;
		
	}
}

// Die Funktion show_map gibt das Array in Form einer 3x3-Matrix aus.
void show_map (void)
{
	for(int x = 0; x < 3; x++)
	{
		for(int y = 0; y < 3; y++)
		{
			printf( "%s ", map [ x ][ y ] );
		}
		printf( "\n" );
	}
}

// Main Method
int main (void)
{
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);
	set_dir(1, 3, N|S|E);
	set_dir(1, 1, N|S|E|W);

	show_map();

	return 0;
}
