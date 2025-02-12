import meshio


def animate_pvds(pvdfile, filename, num_timesteps):

    file1 = open(pvdfile, "w")

    header = '<?xml version="1.0"?>'
    subheader1 = '<VTKFile type="Collection" version="0.1"'
    subheader2 = '         byte_order="LittleEndian"'
    subheader3 = '         compressor="vtkZLibDataCompressor">'

    file1.write(header)
    file1.write("\n")
    file1.write(subheader1)
    file1.write("\n")
    file1.write(subheader2)
    file1.write("\n")
    file1.write(subheader3)
    file1.write("\n")

    start_collection = '  <Collection>'
    file1.write(start_collection)
    file1.write("\n")
    for i in range(num_timesteps):
        timestep = str(i)
        name = f"result/step{i}/step{i}_0.vtu"

        dataset = '    <DataSet timestep="'+timestep+'" group="" part="0"'
        file1.write(dataset)
        file1.write("\n")
        datafile = '             file="'+name+'"/>'
        file1.write(datafile)
        file1.write("\n")

    upfooter = '  </Collection>'
    file1.write(upfooter)
    file1.write("\n")
    footer = '</VTKFile>'
    file1.write(footer)

    file1.close()
