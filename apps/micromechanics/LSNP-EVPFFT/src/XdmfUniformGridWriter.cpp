#include "XdmfUniformGridWriter.h"

XdmfUniformGridWriter::XdmfUniformGridWriter(int N1_, int N2_, int N3_)
  : N1 (N1_)
  , N2 (N2_)
  , N3 (N3_)
  , xdmf (nullptr)
{
}

void XdmfUniformGridWriter::open_file(const char* xdmf_filename)
{
 xdmf = fopen(xdmf_filename, "w");
}

void XdmfUniformGridWriter::write_header(const char* gridName)
{
    fprintf(xdmf, "<?xml version=\"1.0\"?>\n");
    fprintf(xdmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]>\n");
    fprintf(xdmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n");

    fprintf(xdmf, " <Domain>\n");
    fprintf(xdmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n", gridName);
    fprintf(xdmf, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d \"></Topology>\n", N1+1, N2+1, N3+1);
    fprintf(xdmf, "       <Geometry Type=\"ORIGIN_DXDYDZ\">\n");
    fprintf(xdmf, "         <!-- Origin  Z, Y, X -->\n"); // comment
    fprintf(xdmf, "         <DataItem Format=\"XML\" Dimensions=\"3\">%d %d %d</DataItem>\n", 0, 0, 0);
    fprintf(xdmf, "         <!-- DxDyDz (Spacing/Resolution) Z, Y, X -->\n"); // comment
    fprintf(xdmf, "         <DataItem Format=\"XML\" Dimensions=\"3\">%d %d %d</DataItem>\n", 1, 1, 1);
    fprintf(xdmf, "       </Geometry>\n");
}

void XdmfUniformGridWriter::write_attribute(const char* aName, const char* aType,
        const char* aNumberType, const char* aPrecision, const char* aLocation)
{
    fprintf(xdmf, "       <Attribute Name=\"%s\" AttributeType=\"%s\" Center=\"Cell\">\n", aName, aType);
    fprintf(xdmf, "         <DataItem Format=\"HDF\" Dimensions=\"%d %d %d\" NumberType=\"%s\" Precision=\"%s\" >\n", N1, N2, N3, aNumberType, aPrecision);
    fprintf(xdmf, "           %s\n", aLocation);
    fprintf(xdmf, "         </DataItem>\n");
    fprintf(xdmf, "       </Attribute>\n");
}

void XdmfUniformGridWriter::write_footer()
{
    fprintf(xdmf, "   </Grid>\n");
    fprintf(xdmf, " </Domain>\n");
    fprintf(xdmf, "</Xdmf>\n");
}

void XdmfUniformGridWriter::close_file()
{
  if (xdmf != nullptr) {
    fclose(xdmf);
  } 
}

XdmfUniformGridWriter::~XdmfUniformGridWriter()
{
}
