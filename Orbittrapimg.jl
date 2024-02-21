using Images, Interpolations

check_rectangle(a, b, c, d) = (a<b && c<d) || error("[$a $b]x[$c,$d] is not a rectangle")

mutable struct Imageotrap{F, S, T}
  f::F
  img::Union{Matrix{RGB{T}}, Matrix{RGBA{T}}, Matrix{Gray{T}}}
  imgrectgl::NamedTuple{(:a, :b, :c, :d), NTuple{4, S}}  
  fractrectgl::NamedTuple{(:xmin, :xmax, :ymin, :ymax), NTuple{4, S}}
  resfr::Tuple{Int,  Int} #fractal image resolution
end
function Imageotrap(f::F, img::Union{Matrix{RGB{T}}, Matrix{RGBA{T}}, Matrix{Gray{T}}}, 
                    imgrectgl::NamedTuple{(:a, :b, :c, :d), NTuple{4, S}} ,
                    fractrectgl::NamedTuple{(:xmin, :xmax, :ymin, :ymax), NTuple{4, S}}, 
                    resfr::NTuple{2, Int}) where {F<:Function, S<:Real, T<:Real} 
    check_rectangle(imgrectgl...)
    check_rectangle(fractrectgl...)
    return Imageotrap{F, S, T}(f, img, imgrectgl, fractrectgl, resfr) 
end

function get_index(z::Complex, obj::Imageotrap)
    a, b, c, d = obj.imgrectgl
    nrow, ncol = size(obj.img)
    idxrow = floor(Int, 1.5+ (nrow-1)*(d-imag(z))/(d-c)) #added 0.5 for better int approximation
    idxcol = floor(Int, 1.5+ (ncol-1)*(real(z)-a)/(b-a))
    return idxrow, idxcol
end 

function trapped(img, i::Int, j::Int)
       1 <= i <=size(img, 1)  && 1 <= j <= size(img, 2) && (img[i, j] != RGBA(0,0,0,0))     
end 

function iteratef(z::Complex, obj::Imageotrap; maxiter=512) 
    n=0
    i, j = get_index(z, obj)
    while (n < maxiter && abs(z) < 2) && (n < 2 || !trapped(obj.img, i, j ))
        z = obj.f(z)
        n += 1
        i, j = get_index(z, obj;)
        if trapped(obj.img, i, j)
            return obj.img[i, j]
        end
    end
end    

function image2Juliaset(obj::Imageotrap;  maxiter=1024,
                          bgcolor=RGB{N0f8}(1.0, 1.0, 1.0))
    A, B, C, D = obj.fractrectgl
    Nrow, Ncol = obj.resfr
    frimg = fill(bgcolor, 1:Nrow, 1:Ncol)
    for l in axes(frimg, 2)
        for k in  axes(frimg, 1)
            z = Complex(A+(B-A)*(l-1)/(Ncol-1), D-(D-C)*(k-1)/(Nrow-1)) 
            vr = iteratef(z, obj; maxiter=maxiter) 
            if vr != nothing
                frimg[k, l] = vr
            end    
        end
    end
    return interpolate(frimg, BSpline(Linear()), OnGrid())
end;  