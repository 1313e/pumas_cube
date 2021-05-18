// Standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// PUMAS
#include "pumas.h"

// Main function for dumping physics to a binary file
int main(){
    // Obtain name of XML file
    const char *XML_file = "./materials/materials.xml";

    // Obtain length of XML file name
    int len = strlen(XML_file);

    // Obtain filename without extension
    char prefix[len-3];
    strncpy(prefix, &XML_file[0], len-4);
    prefix[len-4] = '\0';

    // Declare physics object
    struct pumas_physics *physics;

    // Create physics
    pumas_physics_create(&physics, PUMAS_PARTICLE_MUON, XML_file,
                         "materials", NULL);

    // Dump to binary file
    char *dump_file;
    sprintf(dump_file, "%s_dump", prefix);
    FILE *file = fopen(dump_file, "wb+");
    pumas_physics_dump(physics, file);
    fclose(file);

    // Release memory
    pumas_physics_destroy(&physics);
}