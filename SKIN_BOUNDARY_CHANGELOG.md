# Winglet skin/boundary changes

Recorded: 2026-10-21

## Problem observed

The first attempt wrote `BPrm->Bndry[4] = SkinT`. This was invalid because
`IritUserMicroTileBndryPrmStruct::Bndry` contains exactly four elements,
indexed 0 through 3. IRIT creates a face skin only when all four corner
thicknesses are positive.

After setting all four values to `0.05`, Winglet geometry generation produced
many negative-Jacobian trivariates while using approximate/global FFD mapping.

## Current repair in `test.c`

1. `PreProcessTile1FaceParam` accepts a four-value array:

   ```c
   const CagdRType SkinThickness[4]
   ```

2. The array is copied into the structure without a loop:

   ```c
   memcpy(BPrm -> Bndry, SkinThickness, sizeof(BPrm -> Bndry));
   ```

3. Lattice type 0 defines uniform arrays:

   ```c
   static const CagdRType NoSkin[4] =
       { 0.0, 0.0, 0.0, 0.0 };
   static const CagdRType UniformSkin[4] =
       { 0.05, 0.05, 0.05, 0.05 };
   ```

4. `UniformSkin` is used only on the body's exterior V faces:

   - V-min when `CBData->TileIdxs[1] == 0`
   - V-max when `CBData->TileIdxs[1] == Body->Counts[1] - 1`

   Interior V faces and all U/W faces use `NoSkin`, avoiding internal solid
   sheets between tiles.

5. `GenerateMicroStructures` currently uses exact local deformation:

   ```c
   const int UseMicroStructFFD = 0;
   ```

   It was previously `1`. Approximate/global FFD mapping generated inverted
   skin trivariates on the strongly curved Winglet deformation map.

## Verification performed

- Reproduced the old geometry failure with the Winglet configuration.
- The old path emitted a large sequence of negative-Jacobian warnings and
  stopped before processing body 2.
- Confirmed the revised `test.c` compiles successfully.
- A freshly rebuilt executable is required before runtime testing.

## If the change must be restored later

Ask to reapply the changes documented in this file. Check that:

- there is no `Bndry[4]` access;
- the helper parameter is `const CagdRType SkinThickness[4]`;
- `NoSkin` and `UniformSkin` are present in lattice type 0;
- only exterior V-min/V-max faces receive `UniformSkin`;
- `UseMicroStructFFD` is `0` for the exact deformation path.
