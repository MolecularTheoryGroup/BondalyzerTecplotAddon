/*
*****************************************************************
*****************************************************************
*******                                                  ********
*******    (C) Copyright 1988-2013 by TECPLOT INC.       *******
*******                                                  ********
*****************************************************************
*****************************************************************
*/


/*
* BondBundlesCreateVolumeZones: For a given BaseZone, cycle through the
* Bonds and create a volume zone contained within the PlsAtomRingCage
* and MnsAtomRingCage surface.
* 1. Find the Bond-PlsAtom and Bond-MnsAtom directions.
* 2. For each point on the Bond-Ring-Cage surface zones, find the
*    intersection of lines projecting in the Bond-MnsAtom and Bond-PlsAtom
*    directions with the MnsAtomRingCage and PlsAtomRingCage surfaces.
* 3. The K=KAve "plane" of the volume zone is the Bond-Ring-Cage surfaces.
*    Successive planes, on either side of the K=KAve plane, are created by
*    evenly spacing points along the projection lines between the Bond-Ring-Cage
*    surface and the intersection points.
*
* Finally, add the new volume zone numbers to the bondBundle structures.
*
* param baseZoneNum
*     Number of the zone for which the topology is computed.
* param volCritPoints
*     Critical point structure for the volume critical points.
* param BondBundlesList
*     ArrayList of pointers to bondBundle data structures for each bond.
*
* return
*     TRUE if successful, FALSE if there were errors.
*/
Boolean_t BondBundlesCreateVolumeZones(const EntIndex_t    baseZoneNum,
                                       const CritPoints_pa volCritPoints,
                                       BondBundles_pa      bondBundles,
                                       double              xMin,
                                       double              xMax,
                                       double              yMin,
                                       double              yMax,
                                       double              zMin,
                                       double              zMax);
