function [] = PlotFmap(C,R,Phi_tilde,phiInd)
%==========================================================================
% Show C, pinv(C) and R
%==========================================================================
figure; 
subplot(1,3,1)
    imagesc(R); colorbar;
    title('${\bf R}$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(1,3,2)
    imagesc(C); colorbar;
    title('${\bf C}$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(1,3,3)    
    imagesc(pinv(C)); colorbar;
    title('${\bf C}^\dagger$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
set(gcf,'Position', [400 400 1500 400])


%==========================================================================
% Test identity
%==========================================================================
PhiTilde_PhiTildeInv = Phi_tilde*pinv(Phi_tilde);
PhiTildeInv_PhiTilde = pinv(Phi_tilde)*Phi_tilde;

figure; 
subplot(1,2,1)
    imagesc(PhiTilde_PhiTildeInv); colorbar;
    title('$\tilde{\Phi} \Phi^\dagger$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(1,2,2)
    imagesc(PhiTildeInv_PhiTilde); colorbar;
    title('$\Phi^\dagger \tilde{\Phi}$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
set(gcf,'Position', [400 400 1000 400])


end