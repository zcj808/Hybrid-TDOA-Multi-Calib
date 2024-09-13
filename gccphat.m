function tau = gccphat(sig, refsig, fs)
    % 计算信号的FFT
    SIG = fft(sig);
    REFSIG = fft(refsig);
    
    % 计算互谱
    R = SIG .* conj(REFSIG);
    
    % 计算GCC-PHAT
    R = R ./ abs(R);
    cc = ifft(R, 'symmetric');
    
    % 找到最大值对应的时间延迟
    [~, max_idx] = max(abs(cc));
    N = length(sig);
    if max_idx > N / 2
        max_idx = max_idx - N;
    end
    tau = max_idx / fs;
end
