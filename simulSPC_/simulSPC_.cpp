// simulSPC_.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include "pch.h"
#include <iostream>
#include "gsl/include/gsl/gsl_integration.h"
#include "Experiment.h"


double f(double x, void * params) {
	double alpha = *(double *)params;
	double f = log(alpha*x) / sqrt(x);
	return f;
}

int main()
{
    std::cout << "Hello World!\n"; 

	gsl_integration_workspace * w
		= gsl_integration_workspace_alloc(1000);

	double result, error;
	double expected = -4.0;
	double alpha = 1.0;

	gsl_function F;
	F.function = &f;
	F.params = &alpha;

	gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000,
		w, &result, &error);

	printf("result          = % .18f\n", result);
	printf("exact result    = % .18f\n", expected);
	printf("estimated error = % .18f\n", error);
	printf("actual error    = % .18f\n", result - expected);
	printf("intervals       = %zu\n", w->size);

	gsl_integration_workspace_free(w);

	return 0;


/*
//TODO read a ini file
Experiment exp;

for(int i=0 ; i < 10000; i++)
	exp.tick();
	
	return 0;
	*/

}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Conseils pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
