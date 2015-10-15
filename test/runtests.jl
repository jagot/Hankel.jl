using Hankel
using Base.Test

using PyCall
pygui(:qt)
using PyPlot

function hankel_test(f1, f2, fun_latex,
                     p, N, R, Va, Na=1000,
                     plot_error=false,
                     fignum=1)
    r = linspace(0,R,Na)
    ν_a = linspace(0,Va,Na)
    figure(fignum,figsize=[12, 17])
    clf()
    for i = 1:length(p)
        subplot(length(p)+1,1,1)
        plot(r, f1(r,p[i]), label = latexstring("\$f_{1,$(p[i])}\$"))
        margins(0.05, 0.1)
        title(latexstring("\$f_{1,p}=$(fun_latex)\$"))
        legend()

        subplot(length(p)+1,1,i+1)
        ax1 = gca()
        ax2 = ax1[:twinx]()
        lines = []
        append!(lines, ax1[:plot](ν_a, abs(f2(ν_a,p[i])),
                                  "k", label="Exact", linewidth = 2.0))
        function test_discrete(trf, name, args...)
            ν_d,f2_d = trf(r -> f1(r,p[i]), p[i], args...)
            append!(lines, ax1[:plot](ν_d, abs(f2_d), "-", label=name))
            if plot_error
                append!(lines, ax2[:semilogy](ν_d, abs(abs(f2(ν_d,p[i]))-abs(f2_d))./abs(maximum(f2_d)), ":",
                                              label="$name error"))
            end
        end

        test_discrete(guizar, "Guizar", R, N)

        ax1[:margins](0.05, 0.1)
        ax1[:set_ylabel]("Function value")
        ax1[:set_title](latexstring("\$f_{2,$(p[i])}\$"))
        ax2[:margins](0.05, 0.1)
        ax2[:set_ylabel]("Dynamic error")
        legend(lines, [l[:get_label]() for l in lines])
    end
    tight_layout()
end

γ = 5
hankel_test((r,p) -> sinc(2π*γ*r),
            function(ν,p)
            function f(ν,p)
            s2 = γ^2 - ν.^2
            if 0<=ν && ν<=γ
            ν^p*cos(p*0.5π)/(2π*γ*√(s2)*(γ*√(s2))^p)
            elseif ν==γ
            Inf
            else
            sin(p*asin(γ./ν))/(2π*γ*√(-s2))
            end
            end
            map(x -> f(x,p), ν)
            end,
            "\\operatorname{sinc}(2\\pi\\gamma r),\\quad \\gamma=$γ",
                0:2, 256, 3.0, 20.0, 1000, true, 1)

a = 0.5
hankel_test((r,p) -> r.^p.*exp(-a*r.^2),
            (ν,p) -> (π/a).^(p+1) * ν.^p .* exp(-π^2*ν.^2/a),
            "r^p\\exp(-ar^2),\\quad a=$a",
            0:4, 256, 20, 7, 1000, true, 2)

hankel_test((r,p) -> ((0 .<= r) .== (r .<= 1))*1.0,
            (ν,p) -> besselj(1,2π*ν)./ν,
            "\\theta(r)",
            0, 128, 20, 15, 1000, true, 3)

# write your own tests here
@test 1 == 1
