function [] = eyediagram(outputpower, time, numb, bitrate, x)

  Tbit = 1e9/bitrate
  time(2) - time(1)
  inter = floor(Tbit/(time(2) - time(1)))
  pause
  i = 1;
  for j = 1 : inter : numel(outputpower)
    if (floor(inter * 1.5) + j < numel(outputpower))
      j
      figure(x)
      plot(time(j : floor(1.5 * inter) + j - 1) - time(j), outputpower(j : floor(1.5 * inter) + j - 1))
      hold on
      i = i + 1
      pause
    end
  end
  figure(x)
  grid on

end
