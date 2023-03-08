using Plots;
using Statistics;

function lin_reg(xdata, ydata)
    x_ = mean(xdata)
    xy_ = mean(xdata .* ydata)
    x2_ = mean(xdata .^ 2)
    x_2 = x_^2
    y_ = mean(ydata)

    beta0 = (x_ * xy_ - x2_ * y_) / (x_2 - x2_)
    beta1 = (x_ * y_ - xy_) / (x_2 - x2_)

    scatter(xdata, ydata, label="Mérési adatok")
    display(plot!(xdata, beta1 .* xdata .+ beta0, label="Illesztett görbe"))

    return beta0, beta1
end


function proba1(xdata, ydata) # f(x) = a*exp(x)+b
    expx = mean(exp.(xdata))
    exp2x = mean(exp.(2 .* xdata))
    expx2 = mean(exp.(xdata)) .^ 2
    y = mean(ydata)
    yexpx = mean(ydata .* exp.(xdata))

    a = (expx .* y - yexpx) / (expx2 - exp2x)
    b = (expx .* yexpx - exp2x .* y) / (expx2 - exp2x)

    scatter(xdata, ydata, label="Mérési adatok")
    display(plot!(xdata, a .* exp.(xdata) .+ b, label="Illesztett görbe"))

    return a, b
end
