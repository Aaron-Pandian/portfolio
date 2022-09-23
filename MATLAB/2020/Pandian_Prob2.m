aapl = load('aapl.txt');
msft = load('msft.txt');
qcom = load('qcom.txt');
wfc = load('wfc.txt');

aaplMean = mean(aapl);
msftMean = mean(msft);
qcomMean = mean(qcom);
wfcMean = mean(wfc);

aapl_msft = ((aapl - aaplMean)'*(msft - msftMean))/ (sqrt((aapl - aaplMean)'*(aapl - aaplMean))* (sqrt((msft - msftMean)'*(msft - msftMean))));
aapl_qcom = ((aapl - aaplMean)'*(qcom - qcomMean))/ (sqrt((aapl - aaplMean)'*(aapl - aaplMean))* (sqrt((qcom - qcomMean)'*(qcom - qcomMean))));
aapl_wfc = ((aapl - aaplMean)'*(wfc - wfcMean))/ (sqrt((aapl - aaplMean)'*(aapl - aaplMean))* (sqrt((wfc - wfcMean)'*(wfc - wfcMean))));

disp(aapl_msft);
disp(aapl_qcom);
disp(aapl_wfc);

hold on;
plot(aapl);
plot(msft);
plot(qcom);
plot(wfc);
