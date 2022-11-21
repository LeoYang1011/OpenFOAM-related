IOobjectList objects(mesh, runTime.timeName());
const wordList& fieldsName = spanwise.fields();

List<wordRe> fieldRes(fieldsName.size());
for (label i=0;i<fieldsName.size();i++)
{
    wordRe fieldRe(fieldsName[i]);
    fieldRes[i] = fieldRe;
}
wordRes fRes(fieldRes);
IOobjectList fobjects = objects.lookup(fRes);
  
PtrList<volScalarField> vsFlds;
ReadFields(mesh, fobjects, vsFlds);

PtrList<volVectorField> vvFlds;
ReadFields(mesh, fobjects, vvFlds);

PtrList<volSphericalTensorField> vstFlds;
ReadFields(mesh, fobjects, vstFlds);

PtrList<volSymmTensorField> vsymtFlds;
ReadFields(mesh, fobjects, vsymtFlds);

PtrList<volTensorField> vtFlds;
ReadFields(mesh, fobjects, vtFlds);
