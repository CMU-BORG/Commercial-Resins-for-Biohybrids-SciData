function [mu_SF, sig_SF] = RoundToSigFigs(mu,sig)
%{
    This function returns the input mean (mu) and standard deviation (sig)
    to the appropriate significant significant figures. Here, the number
    place to round to is the first non-zero digit in the error, unless that
    first digit is 1, where we round to one more digit.
%}

sig_str = num2str(sig);
mu_str = num2str(mu);

dp_sig = find(sig_str=='.');
round_dig = find((sig_str ~= '0')&(sig_str ~='.')&(sig_str ~='-'),1);
if length(dp_sig)==0
    dp_sig = length(sig_str)+1;
end

dp_mu = find(mu_str=='.');
round_dig_mu = find((mu_str ~= '0')&(mu_str ~='.'),1);
if length(dp_mu)==0
    dp_mu = length(mu_str)+1;
end

if sig_str(round_dig)=='1'
    round_dig = round_dig+1;
end


from_dp = round_dig - dp_sig;
from_dp_mu = round_dig_mu - dp_mu;

if from_dp<0
    from_dp = from_dp + 1;
end
if from_dp_mu<0
    from_dp_mu = from_dp_mu + 1;
end

sig_SF = round((10^from_dp)*sig) / (10^from_dp);

if from_dp_mu >= from_dp
    mu_SF = round((10^from_dp_mu)*mu) / (10^from_dp_mu);
else
    mu_SF = round((10^from_dp)*mu) / (10^from_dp);
end

% fprintf("%s +/- %s\n",num2str(mu_SF),num2str(sig_SF))

end