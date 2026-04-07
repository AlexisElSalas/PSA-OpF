#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvm.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "pimpleControl.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"
#include "OSstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
	//#include "basicMultiComponentMixture.H"

    while (runTime.loop())
    {
		
		// ALMACENAR TIEMPO ANTERIOR PARA VARIABLES ALGEBRAICAS
       // qCO2.storeOldTime();
       // qN2.storeOldTime();
		
        while (pimple.loop())
        {
            // 1. Cinética de Adsorción y términos fuente (Eq. 24 de Henry)
            #include "adsorptionEqn.H"
            
            // 2. Conservación de Especies (Eq. 19 de Henry)
            #include "YEqn.H"
            
            // NOTA: Se eliminó #include "TsEqn.H" para cumplir con el modelo LTE
            
            // 3. Equilibrio Térmico Local - Gas y Sólido (Eq. 22 de Henry)
            #include "TEqn.H"
            
            // 4. Predictor de Momento
            #include "UEqn.H"
            
            // 5. Corrector de Presión / Continuidad
            #include "pEqn.H"
        }

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}