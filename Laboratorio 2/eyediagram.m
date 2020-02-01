function [] = eyediagram(outputpower, bitrate, time, x)

  T_bit = 1/bitrate;
  numb = ceil(numel(time)/T_bit);

  for j = 0 : numb : numel(time)
    figure(99)
    plot(time(j : numb + j), outputpower(j : numb + j))
    hold on
  end
end
