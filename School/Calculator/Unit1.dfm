object Form1: TForm1
  Left = 889
  Top = 174
  BorderIcons = [biSystemMenu, biMinimize]
  BorderStyle = bsSingle
  Caption = #1050#1072#1083#1100#1082#1091#1083#1103#1090#1086#1088
  ClientHeight = 222
  ClientWidth = 198
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object btn1: TButton
    Left = 8
    Top = 32
    Width = 33
    Height = 33
    Caption = '1'
    TabOrder = 0
    OnClick = btn1Click
  end
  object Panel1: TPanel
    Left = 0
    Top = 0
    Width = 198
    Height = 25
    Align = alTop
    BevelOuter = bvNone
    TabOrder = 1
    object edtInput: TEdit
      Left = 0
      Top = 0
      Width = 201
      Height = 21
      BiDiMode = bdLeftToRight
      ParentBiDiMode = False
      ReadOnly = True
      TabOrder = 0
      Text = '0'
    end
  end
  object btn2: TButton
    Left = 56
    Top = 32
    Width = 33
    Height = 33
    Caption = '2'
    TabOrder = 2
    OnClick = btn1Click
  end
  object btn3: TButton
    Left = 104
    Top = 32
    Width = 33
    Height = 33
    Caption = '3'
    TabOrder = 3
    OnClick = btn1Click
  end
  object btn4: TButton
    Left = 8
    Top = 80
    Width = 33
    Height = 33
    Caption = '4'
    TabOrder = 4
    OnClick = btn1Click
  end
  object btn5: TButton
    Left = 56
    Top = 80
    Width = 33
    Height = 33
    Caption = '5'
    TabOrder = 5
    OnClick = btn1Click
  end
  object btn6: TButton
    Left = 104
    Top = 80
    Width = 33
    Height = 33
    Caption = '6'
    TabOrder = 6
    OnClick = btn1Click
  end
  object btn7: TButton
    Left = 8
    Top = 128
    Width = 33
    Height = 33
    Caption = '7'
    TabOrder = 7
    OnClick = btn1Click
  end
  object btn8: TButton
    Left = 56
    Top = 128
    Width = 33
    Height = 33
    Caption = '8'
    TabOrder = 8
    OnClick = btn1Click
  end
  object btn9: TButton
    Left = 104
    Top = 128
    Width = 33
    Height = 33
    Caption = '9'
    TabOrder = 9
    OnClick = btn1Click
  end
  object btn_t: TButton
    Left = 8
    Top = 176
    Width = 33
    Height = 33
    Caption = '.'
    TabOrder = 10
  end
  object btn0: TButton
    Left = 56
    Top = 176
    Width = 33
    Height = 33
    Caption = '0'
    TabOrder = 11
    OnClick = btn1Click
  end
  object btn_r: TButton
    Left = 104
    Top = 176
    Width = 33
    Height = 33
    Caption = '='
    TabOrder = 12
    OnClick = btn_rClick
  end
  object btn_p: TButton
    Left = 152
    Top = 32
    Width = 33
    Height = 33
    Caption = '+'
    TabOrder = 13
    OnClick = btn_pClick
  end
  object btn_m: TButton
    Left = 152
    Top = 80
    Width = 33
    Height = 33
    Caption = '-'
    TabOrder = 14
    OnClick = btn_pClick
  end
  object btn_d: TButton
    Left = 152
    Top = 128
    Width = 33
    Height = 33
    Caption = '/'
    TabOrder = 15
    OnClick = btn_pClick
  end
  object btn_u: TButton
    Left = 152
    Top = 176
    Width = 33
    Height = 33
    Caption = '*'
    TabOrder = 16
    OnClick = btn_pClick
  end
end
