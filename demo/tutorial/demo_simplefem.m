%% Solving stochastic PDEs
% This demonstration show how to solve a parametric stochastic PDE using
% the GPC methods in sglib. 'Parametric' here means, that the uncertainties
% in PDE are given by a finite set of parameters, and not by a stochastic
% field described by infinitely many stochastic variable, which has to be
% (stochastically) discretised beforehand.

publishing_defaults

%% The deterministic problem
% Lets look at the following boundary value problem 
%
% $$-\frac{\partial}{\partial x}\left(a(x) \frac{\partial}{\partial x}\right)=1 \textrm{ on } 
% \mathcal{D}=[0,1]$$
% 
% where a(x) is a coefficient field given by 
%
% $$a(x) = \left\{\!\! \begin{array}{ll} a_1, & x<0.5 \\ a_2, & x\geq0.5 \end{array} \right. $$
% 

%%
% The code for solving this problem on a grid with N=101 nodes has been
% put into |diffusion_1d_solve| which takes two parameters. Here is an
% example solving the problem for a_1=2 and a_2=100.

a1 = 2; a2 = 100;
[u, a, pos]=diffusion_1d_complete_solve(a1, a2);
plot(pos, u); 





%% Stochastic description of the parameters
% Now suppose the parameters a_1 and a_2 are uncertain and described by
% by random variables. Suppose a_1 follows a Beta[1.2,2] distribution,
% which is shift and scaled such that it has support [0.5, 5]. In sglib
% this can be easily specified within the gendist frame work, using the
% |gendist_create| and |gendist_fix_bounds| functions (see
% demo_distributions). 

a1_dist = gendist_create('beta', {1.2, 2});
a1_dist = gendist_fix_bounds(a1_dist, 0.5, 5);

%%
% The probability distribution function (pdf) and the cumulative
% distribution function (cdf) can be computed with |gendist_cdf| and
% |gendist_pdf| function (ok, it's now all in the plot function)

density_plot(a1_dist, 'type', 'both')
legend( 'pdf a1', 'cdf a1');

%%
% The second parameter a_2 follows a shifted and scaled Beta[0.3, 0.6]
% distribution with support on [50, 150]. The distribution of a_2 has
% sharp singularities at both ends.

a2_dist = gendist_create('beta', {0.6, 0.3});
a2_dist = gendist_fix_bounds(a2_dist, 50, 150);

density_plot(a2_dist, 'type', 'both')
legend( 'pdf a2', 'cdf a2');



%% GPC approximation of the parameters
% For some methods we could directly use the distributions given above but
% here we want to use their gpc approximation.

%%
% For $a_1$ the distribution looks quite a bit like a (Wigner) semicircle
% distribution for which we can use the combination Semicircle/ChebyshevU,
% specified in sglib by the letter 'u' (small 'u' for the normalised
% polynomials). With the option 'varerr' we instruct |gpc_param_expand| to
% choose such a degree of the gpc polynomial, that the error in variance is
% below 0.001. The option 'fixvar' rescale the gpc coefficients such that
% the variance is fully matched (possibly incurring a worse match in higher
% moments).
% 

[a1_alpha, V1, err] = gpc_param_expand(a1_dist, 'u', 'varerr', 0.001, 'fixvar', true);

%%
% Compare the first four moments of the true distribution of a1 and its GPC
% approximation

[mean,var,skew,kurt]=gendist_moments(a1_dist);
fprintf('Moments (a1,true):\nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)
[mean,var,skew,kurt]=gpc_moments(a1_alpha, V1);
fprintf('Moments (a1,gpc): \nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)

%%
% Draw some samples from the underlying distribution, then plot a kernel
% density estimate of the gpc approximation of a1 and compare to the true
% distribution.
a1_samples = gpc_evaluate(a1_alpha, V1, gpcgerm_sample(V1, 100000)); 

density_plot(a1_samples, 'type', 'kernel', 'rug', true, 'max_rug', inf);
density_plot(a1_dist, 'hold', true);
legend('gpc approx. (kde)', 'samples', 'exact density');

%%
% It can be seen in the plot, that ... (somehow a bit contradictory with
% the kurtosis...hmmm)

%% 
% For the $a_2$ the distribution looks quite a bit like an arcsine
% distribution for which we can use the combination Arcsine/ChebyshevT,
% specified in sglib by the letter 't'.

[a2_alpha, V2] = gpc_param_expand(a2_dist, 't', 'varerr', 0.001, 'fixvar', true);


%%
% Compare the first four moments of the true distribution of a1 and its GPC
% approximation
[mean,var,skew,kurt]=gendist_moments(a2_dist);
fprintf('Moments (a2,true):\nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)
[mean,var,skew,kurt]=gpc_moments(a2_alpha, V2);
fprintf('Moments (a2,gpc): \nmean=%g, var=%g, skew=%g, kurt=%g\n', mean, var, skew, kurt)

%%
% Plot a kernel density estimate of the gpc approximation of a2 and compare
% to the true distribution
a2_samples = gpc_evaluate(a2_alpha, V2, gpcgerm_sample(V2, 100000));

density_plot(a2_samples, 'type', 'kernel', 'kde_sig', 0.025, 'rug', true, 'max_rug', 300);
density_plot(a2_dist, 'hold', true);
legend('gpc approx. (kde)', 'samples', 'exact density');




%% Combining the GPC spaces
% The gpc spaces V1 anv V2 constructed for the parameters a1 and a2 are
% currently disjoined spaces. To use them in any of the further algorithms
% (except maybe Monte Carlo), we need to combine them into larger space,
% i.e. take the product space). This can be achieved by the
% |gpc_combine_inputs| function (which can be used much more generally than
% the name implies.). This will not only form the product space, but also
% map the coefficients to the right places in the new combined coefficient
% field.

[a_i_alpha, V_a] = gpc_combine_inputs(a1_alpha, V1, a2_alpha, V2);

%%
% If we sample from the germ distribution, we can see that it's indeed now
% two-dimensional and has a Semicircle x Arcsine distribution
xi = gpcgerm_sample(V_a, 100000);
plot(xi(1,:), xi(2,:), '.', 'MarkerSize', 1); axis equal;

%%
% Sampling from the parameters (now with quasi Monte Carlo) gives the
% following (note that the arcsine distribution goes vertically with the
% larger spike on top, and the semicircle distribtion goes horizontally
% with the bump more to the left).
a_i_samples = gpc_sample(a_i_alpha, V_a, 100000, 'mode', 'qmc');
plot(a_i_samples(1,:), a_i_samples(2,:), '.', 'MarkerSize', 1);

%% Monte Carlo
% Now, as a first step, we can do a Monte-Carlo simulation of our model.
% Instead of making statistics, we'll just plot a bunch of samples.

N = 30;
a_i_samples = gpc_sample(a_i_alpha, V_a, N, 'mode', 'qmc');
u_samples = zeros(N, size(pos,2));
[state, info] = diffusion_1d_init();
for i=1:N
    u=diffusion_1d_solve(state, a_i_samples(:,i));
    u_samples(i,:) = u;
end
plot(pos, u_samples)

%% And now for Galerkin...
V_u = gpcbasis_create(V_a, 'p', 5);

A_i = gpc_multiplication_matrices(a_i_alpha, V_a, V_u);
subplot(1,2,1); spy(A_i{1})
subplot(1,2,2); spy(A_i{2})

