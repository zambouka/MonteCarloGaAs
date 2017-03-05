#include <iostream>
#include "electron.h"
#include "mainArea.h"
#include <fstream>

using namespace std;
//test impactionization
int main()
{
    ofstream dV("driftVelvsField.txt");
    Electron elctron[5000];
    for (size_t i = 0; i < 1; i++)
    {
        mainArea a(350000);

        for (size_t i = 0; i < 5000; i++)
        {
            a.init(elctron[i]);
        }

        for (size_t j = 0; j < 20000; j++)
        {
            for (size_t i = 0; i < 5000; i++)
            {
                a.process(elctron[i]);
            }
            if (j>4000)
            {
                for (size_t i = 0; i < 5000; i++)
                {
                    a.analysis(elctron[i]);
                }
                a.write();
                double v = a.driftVY(5000);
                dV << v << endl;
                a.clearOut();
            }
        }
    }
    return 0;
}