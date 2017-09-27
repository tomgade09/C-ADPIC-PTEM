def readFiles(filename1, filename2, filename3, filename4, filename5, filename6):
    fn1 = ctypes.create_string_buffer(bytes(filename1, encoding='utf-8'))
    fn2 = ctypes.create_string_buffer(bytes(filename2, encoding='utf-8'))
    fn3 = ctypes.create_string_buffer(bytes(filename3, encoding='utf-8'))
    fn4 = ctypes.create_string_buffer(bytes(filename4, encoding='utf-8'))
    fn5 = ctypes.create_string_buffer(bytes(filename5, encoding='utf-8'))
    fn6 = ctypes.create_string_buffer(bytes(filename6, encoding='utf-8'))
    e_vpara_c = cppDLL.readDblBin(fn1, 100352)
    e_vperp_c = cppDLL.readDblBin(fn2, 100352)
    e_z_c = cppDLL.readDblBin(fn3, 100352)
    i_vpara_c = cppDLL.readDblBin(fn4, 100352)
    i_vperp_c = cppDLL.readDblBin(fn5, 100352)
    i_z_c = cppDLL.readDblBin(fn6, 100352)

    e_vpara = []
    e_vperp = []
    e_z = []
    i_vpara = []
    i_vperp = []
    i_z = []

    for i in range(100352):
        if ((e_z_c[i] < 10) and (e_z_c[i] > (2.0e6 + 6.371e6) / 6.371e6)):
            e_vpara.append(e_vpara_c[i])
            e_vperp.append(e_vperp_c[i])
            e_z.append(e_z_c[i])
        
        if ((i_z_c[i] < 10) and (i_z_c[i] > (2.0e6 + 6.371e6) / 6.371e6)):
            i_vpara.append(i_vpara_c[i])
            i_vperp.append(i_vperp_c[i])
            i_z.append(i_z_c[i])

    return e_vpara, e_vperp, e_z, i_vpara, i_vperp, i_z, [0.], [0.], [0.]