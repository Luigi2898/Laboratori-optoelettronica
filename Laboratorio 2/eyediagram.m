function [] = eyediagram(outputpower, time, numb, bitrate, x)

  Tbit = 1e9/bitrate;
  inter = ceil(Tbit/(time(2) - time(1)));
  %pause
  i = 1;
  for j = 1 : inter : numel(outputpower)
    if (inter + j < numel(outputpower))
      figure(x)
      plot([1 : 1 : floor(1.5 * inter)], outputpower(j : floor(1.5 * inter) + j - 1))
      hold on
      i = i + 1
    end
  end
  figure(x)
  grid on

end

%FIXME partire dalla transizione forse, non da dove inizia
