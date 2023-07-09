@compile_workload begin
    m0 = Multicomplex(1.)
    m1 = Multicomplex(m0, m0)
    m2 = Multicomplex(m1, m1)
    m3 = Multicomplex{3}(SVector{8}(float(1:8)))
    m4 = Multicomplex{4}(SVector{16}(float(1:16)))

    x1 = m0 + m1 + m2 + m3 + m4
    x2 = 2*m0
    x3 = m0 * m1
    x4 = x3 * m2
    x5 = x4 * m3
    x6 = x5 * m4
    x7 = exp(m0)
    x8 = exp(m1)
    x9 = exp(m2)
    x10= exp(m3)
    # x11= exp(m4)
    x12 = m1 / m0
    x13 = m2 / m1
    x14 = m3 / m2
    x15 = m4 / m3

    io = IOBuffer()
    show(io, m0)
    show(io, m1)
    show(io, m2)
    show(io, m3)
    show(io, m4)
end