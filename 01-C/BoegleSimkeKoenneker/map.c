#include <stdio.h>

#define HEIGHT 3
#define WIDTH 3

// Definieren Sie ein enum cardd
typedef enum {
	N=1<<0,
	E=1<<1,
	S=1<<2,
	W=1<<3
} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[HEIGHT][WIDTH];

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
	if (x >= 0 && x < HEIGHT &&
		y >= 0 && y < WIDTH &&
		dir >= 0 && dir <= (N|E|S|W) &&
		!(dir&N && dir&S) &&
		!(dir&E && dir&W))
	{
		map[x][y] = dir;
	}
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
	for (int x = 0; x < HEIGHT; x++)
	{
		for (int y = 0; y < WIDTH; y++)
		{
			char *dir;
			switch (map[x][y])
			{
			case N:
				dir = "N";
				break;
			case E:
				dir = "E";
				break;
			case S:
				dir = "S";
				break;
			case W:
				dir = "W";
				break;
			case N|E:
				dir = "NE";
				break;
			case N|W:
				dir = "NW";
				break;
			case S|E:
				dir = "SE";
				break;
			case S|W:
				dir = "SW";
				break;
			default:
				dir = "0";
			}

			if (y == 0)
			{
				printf("%-2s ", dir);
			}
			else if (y == WIDTH-1)
			{
				printf(" % 2s\n", dir);
			}
			else
			{
				printf(" %-2s", dir);
			}
		}
	}
}

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
