% ================================================================= %
% 3差分解析スクリプト（修正版：複素数・レンジ最適化）
% ================================================================= %

% 1. ファイル選択
[file_led, path_led] = uigetfile('*.mat', '【1/2】LEDあり（80Hz等）の.matを選択');
if isequal(file_led, 0), return; end
[file_noled, path_noled] = uigetfile('*.mat', '【2/2】LEDなし（0Hz等）の.matを選択', path_led);
if isequal(file_noled, 0), return; end

% 2. データのロード
L = load(fullfile(path_led, file_led));
N = load(fullfile(path_noled, file_noled));

% --- 差分の計算（ここが修正ポイント） ---
% ERSPの差
diff_ersp = L.ersp - N.ersp;

% ITCは複素数なので、まず「絶対値（強さ）」にしてから引き算する
diff_itc  = abs(L.itc) - abs(N.itc); 

% ERS（絶対パワー）の差
diff_ers   = L.ers - N.ers;

% --- 描画設定 ---
fig = figure('Color', 'w', 'Position', [100, 50, 600, 900]); 

targets = {diff_ersp, diff_itc, diff_ers};
titles  = {'ERSP Difference (dB)', 'ITC Difference (Magnitude)', 'ERS Difference (dB)'};

% スケール設定（白紙にならないよう、データの最大値に合わせて自動調整する機能を追加）
% ITCは 0.3, 他は 1.5 をデフォルトにしていますが、データが小さければ自動で縮小します
limits  = [max(abs(diff_ersp(:)))*1.0, max(abs(diff_itc(:)))*0.8, max(abs(diff_ers(:)))*3];

for i = 1:3
    subplot(3, 1, i);
    imagesc(L.times, L.freqs, targets{i});
    ylim([0 200])
    set(gca, 'YDir', 'normal', 'FontSize', 10, 'FontWeight', 'bold');
    
    colormap(gca, jet);
    cb = colorbar;
    
    % 差分が極端に小さい場合、白紙に見えるのを防ぐため最小レンジを確保
    current_lim = max(limits(i), 0.05); 
    caxis([-current_lim, current_lim]); 
    
    title(titles{i}, 'FontSize', 12);
    ylabel('Frequency (Hz)');
    
    if i == 3, xlabel('Time (ms)'); else, set(gca, 'XTickLabel', []); end
end

[~, name_led] = fileparts(file_led);
sgtitle(['Difference Analysis: ', name_led], 'Interpreter', 'none', 'FontWeight', 'bold');