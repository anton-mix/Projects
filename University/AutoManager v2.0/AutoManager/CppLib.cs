using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.Runtime.InteropServices;


namespace AutoManager {
	enum WindowsMessages : int {
		 WM_ACTIVATE = 0x0006,
		 WM_ACTIVATEAPP = 0x001C,
		 WM_AFXFIRST = 0x0360,
		 WM_AFXLAST = 0x037F,
		 WM_APP = 0x8000,
		 WM_ASKCBFORMATNAME = 0x030C,
		 WM_CANCELJOURNAL = 0x004B,
		 WM_CANCELMODE = 0x001F,
		 WM_CAPTURECHANGED = 0x0215,
		 WM_CHANGECBCHAIN = 0x030D,
		 WM_CHANGEUISTATE = 0x0127,
		 WM_CHAR = 0x0102,
		 WM_CHARTOITEM = 0x002F,
		 WM_CHILDACTIVATE = 0x0022,
		 WM_CLEAR = 0x0303,
		 WM_CLOSE = 0x0010,
		 WM_COMMAND = 0x0111,
		 WM_COMPACTING = 0x0041,
		 WM_COMPAREITEM = 0x0039,
		 WM_CONTEXTMENU = 0x007B,
		 WM_COPY = 0x0301,
		 WM_COPYDATA = 0x004A,
		 WM_CREATE = 0x0001,
		 WM_CTLCOLORBTN = 0x0135,
		 WM_CTLCOLORDLG = 0x0136,
		 WM_CTLCOLOREDIT = 0x0133,
		 WM_CTLCOLORLISTBOX = 0x0134,
		 WM_CTLCOLORMSGBOX = 0x0132,
		 WM_CTLCOLORSCROLLBAR = 0x0137,
		 WM_CTLCOLORSTATIC = 0x0138,
		 WM_CUT = 0x0300,
		 WM_DEADCHAR = 0x0103,
		 WM_DELETEITEM = 0x002D,
		 WM_DESTROY = 0x0002,
		 WM_DESTROYCLIPBOARD = 0x0307,
		 WM_DEVICECHANGE = 0x0219,
		 WM_DEVMODECHANGE = 0x001B,
		 WM_DISPLAYCHANGE = 0x007E,
		 WM_DRAWCLIPBOARD = 0x0308,
		 WM_DRAWITEM = 0x002B,
		 WM_DROPFILES = 0x0233,
		 WM_ENABLE = 0x000A,
		 WM_ENDSESSION = 0x0016,
		 WM_ENTERIDLE = 0x0121,
		 WM_ENTERMENULOOP = 0x0211,
		 WM_ENTERSIZEMOVE = 0x0231,
		 WM_ERASEBKGND = 0x0014,
		 WM_EXITMENULOOP = 0x0212,
		 WM_EXITSIZEMOVE = 0x0232,
		 WM_FONTCHANGE = 0x001D,
		 WM_GETDLGCODE = 0x0087,
		 WM_GETFONT = 0x0031,
		 WM_GETHOTKEY = 0x0033,
		 WM_GETICON = 0x007F,
		 WM_GETMINMAXINFO = 0x0024,
		 WM_GETOBJECT = 0x003D,
		 WM_GETTEXT = 0x000D,
		 WM_GETTEXTLENGTH = 0x000E,
		 WM_HANDHELDFIRST = 0x0358,
		 WM_HANDHELDLAST = 0x035F,
		 WM_HELP = 0x0053,
		 WM_HOTKEY = 0x0312,
		 WM_HSCROLL = 0x0114,
		 WM_HSCROLLCLIPBOARD = 0x030E,
		 WM_ICONERASEBKGND = 0x0027,
		 WM_IME_CHAR = 0x0286,
		 WM_IME_COMPOSITION = 0x010F,
		 WM_IME_COMPOSITIONFULL = 0x0284,
		 WM_IME_CONTROL = 0x0283,
		 WM_IME_ENDCOMPOSITION = 0x010E,
		 WM_IME_KEYDOWN = 0x0290,
		 WM_IME_KEYLAST = 0x010F,
		 WM_IME_KEYUP = 0x0291,
		 WM_IME_NOTIFY = 0x0282,
		 WM_IME_REQUEST = 0x0288,
		 WM_IME_SELECT = 0x0285,
		 WM_IME_SETCONTEXT = 0x0281,
		 WM_IME_STARTCOMPOSITION = 0x010D,
		 WM_INITDIALOG = 0x0110,
		 WM_INITMENU = 0x0116,
		 WM_INITMENUPOPUP = 0x0117,
		 WM_INPUTLANGCHANGE = 0x0051,
		 WM_INPUTLANGCHANGEREQUEST = 0x0050,
		 WM_KEYDOWN = 0x0100,
		 WM_KEYFIRST = 0x0100,
		 WM_KEYLAST = 0x0108,
		 WM_KEYUP = 0x0101,
		 WM_KILLFOCUS = 0x0008,
		 WM_LBUTTONDBLCLK = 0x0203,
		 WM_LBUTTONDOWN = 0x0201,
		 WM_LBUTTONUP = 0x0202,
		 WM_MBUTTONDBLCLK = 0x0209,
		 WM_MBUTTONDOWN = 0x0207,
		 WM_MBUTTONUP = 0x0208,
		 WM_MDIACTIVATE = 0x0222,
		 WM_MDICASCADE = 0x0227,
		 WM_MDICREATE = 0x0220,
		 WM_MDIDESTROY = 0x0221,
		 WM_MDIGETACTIVE = 0x0229,
		 WM_MDIICONARRANGE = 0x0228,
		 WM_MDIMAXIMIZE = 0x0225,
		 WM_MDINEXT = 0x0224,
		 WM_MDIREFRESHMENU = 0x0234,
		 WM_MDIRESTORE = 0x0223,
		 WM_MDISETMENU = 0x0230,
		 WM_MDITILE = 0x0226,
		 WM_MEASUREITEM = 0x002C,
		 WM_MENUCHAR = 0x0120,
		 WM_MENUCOMMAND = 0x0126,
		 WM_MENUDRAG = 0x0123,
		 WM_MENUGETOBJECT = 0x0124,
		 WM_MENURBUTTONUP = 0x0122,
		 WM_MENUSELECT = 0x011F,
		 WM_MOUSEACTIVATE = 0x0021,
		 WM_MOUSEFIRST = 0x0200,
		 WM_MOUSEHOVER = 0x02A1,
		 WM_MOUSELAST = 0x020D,
		 WM_MOUSELEAVE = 0x02A3,
		 WM_MOUSEMOVE = 0x0200,
		 WM_MOUSEWHEEL = 0x020A,
		 WM_MOUSEHWHEEL = 0x020E,
		 WM_MOVE = 0x0003,
		 WM_MOVING = 0x0216,
		 WM_NCACTIVATE = 0x0086,
		 WM_NCCALCSIZE = 0x0083,
		 WM_NCCREATE = 0x0081,
		 WM_NCDESTROY = 0x0082,
		 WM_NCHITTEST = 0x0084,
		 WM_NCLBUTTONDBLCLK = 0x00A3,
		 WM_NCLBUTTONDOWN = 0x00A1,
		 WM_NCLBUTTONUP = 0x00A2,
		 WM_NCMBUTTONDBLCLK = 0x00A9,
		 WM_NCMBUTTONDOWN = 0x00A7,
		 WM_NCMBUTTONUP = 0x00A8,
		 WM_NCMOUSEHOVER = 0x02A0,
		 WM_NCMOUSELEAVE = 0x02A2,
		 WM_NCMOUSEMOVE = 0x00A0,
		 WM_NCPAINT = 0x0085,
		 WM_NCRBUTTONDBLCLK = 0x00A6,
		 WM_NCRBUTTONDOWN = 0x00A4,
		 WM_NCRBUTTONUP = 0x00A5,
		 WM_NCXBUTTONDBLCLK = 0x00AD,
		 WM_NCXBUTTONDOWN = 0x00AB,
		 WM_NCXBUTTONUP = 0x00AC,
		 WM_NCUAHDRAWCAPTION = 0x00AE,
		 WM_NCUAHDRAWFRAME = 0x00AF,
		 WM_NEXTDLGCTL = 0x0028,
		 WM_NEXTMENU = 0x0213,
		 WM_NOTIFY = 0x004E,
		 WM_NOTIFYFORMAT = 0x0055,
		 WM_NULL = 0x0000,
		 WM_PAINT = 0x000F,
		 WM_PAINTCLIPBOARD = 0x0309,
		 WM_PAINTICON = 0x0026,
		 WM_PALETTECHANGED = 0x0311,
		 WM_PALETTEISCHANGING = 0x0310,
		 WM_PARENTNOTIFY = 0x0210,
		 WM_PASTE = 0x0302,
		 WM_PENWINFIRST = 0x0380,
		 WM_PENWINLAST = 0x038F,
		 WM_POWER = 0x0048,
		 WM_POWERBROADCAST = 0x0218,
		 WM_PRINT = 0x0317,
		 WM_PRINTCLIENT = 0x0318,
		 WM_QUERYDRAGICON = 0x0037,
		 WM_QUERYENDSESSION = 0x0011,
		 WM_QUERYNEWPALETTE = 0x030F,
		 WM_QUERYOPEN = 0x0013,
		 WM_QUEUESYNC = 0x0023,
		 WM_QUIT = 0x0012,
		 WM_RBUTTONDBLCLK = 0x0206,
		 WM_RBUTTONDOWN = 0x0204,
		 WM_RBUTTONUP = 0x0205,
		 WM_RENDERALLFORMATS = 0x0306,
		 WM_RENDERFORMAT = 0x0305,
		 WM_SETCURSOR = 0x0020,
		 WM_SETFOCUS = 0x0007,
		 WM_SETFONT = 0x0030,
		 WM_SETHOTKEY = 0x0032,
		 WM_SETICON = 0x0080,
		 WM_SETREDRAW = 0x000B,
		 WM_SETTEXT = 0x000C,
		 WM_SETTINGCHANGE = 0x001A,
		 WM_SHOWWINDOW = 0x0018,
		 WM_SIZE = 0x0005,
		 WM_SIZECLIPBOARD = 0x030B,
		 WM_SIZING = 0x0214,
		 WM_SPOOLERSTATUS = 0x002A,
		 WM_STYLECHANGED = 0x007D,
		 WM_STYLECHANGING = 0x007C,
		 WM_SYNCPAINT = 0x0088,
		 WM_SYSCHAR = 0x0106,
		 WM_SYSCOLORCHANGE = 0x0015,
		 WM_SYSCOMMAND = 0x0112,
		 WM_SYSDEADCHAR = 0x0107,
		 WM_SYSKEYDOWN = 0x0104,
		 WM_SYSKEYUP = 0x0105,
		 WM_TCARD = 0x0052,
		 WM_TIMECHANGE = 0x001E,
		 WM_TIMER = 0x0113,
		 WM_UNDO = 0x0304,
		 WM_UNINITMENUPOPUP = 0x0125,
		 WM_USER = 0x0400,
		 WM_USERCHANGED = 0x0054,
		 WM_VKEYTOITEM = 0x002E,
		 WM_VSCROLL = 0x0115,
		 WM_VSCROLLCLIPBOARD = 0x030A,
		 WM_WINDOWPOSCHANGED = 0x0047,
		 WM_WINDOWPOSCHANGING = 0x0046,
		 WM_WININICHANGE = 0x001A,
		 WM_XBUTTONDBLCLK = 0x020D,
		 WM_XBUTTONDOWN = 0x020B,
		 WM_XBUTTONUP = 0x020C
}

	enum VirtualKeyStates : int {
		VK_LBUTTON = 0x01,
		VK_RBUTTON = 0x02,
		VK_CANCEL = 0x03,
		VK_MBUTTON = 0x04,
		//
		VK_XBUTTON1 = 0x05,
		VK_XBUTTON2 = 0x06,
		//
		VK_BACK = 0x08,
		VK_TAB = 0x09,
		//
		VK_CLEAR = 0x0C,
		VK_RETURN = 0x0D,
		//
		VK_SHIFT = 0x10,
		VK_CONTROL = 0x11,
		VK_MENU = 0x12,
		VK_PAUSE = 0x13,
		VK_CAPITAL = 0x14,
		//
		VK_KANA = 0x15,
		VK_HANGEUL = 0x15,  /* old name - should be here for compatibility */
		VK_HANGUL = 0x15,
		VK_JUNJA = 0x17,
		VK_FINAL = 0x18,
		VK_HANJA = 0x19,
		VK_KANJI = 0x19,
		//
		VK_ESCAPE = 0x1B,
		//
		VK_CONVERT = 0x1C,
		VK_NONCONVERT = 0x1D,
		VK_ACCEPT = 0x1E,
		VK_MODECHANGE = 0x1F,
		//
		VK_SPACE = 0x20,
		VK_PRIOR = 0x21,
		VK_NEXT = 0x22,
		VK_END = 0x23,
		VK_HOME = 0x24,
		VK_LEFT = 0x25,
		VK_UP = 0x26,
		VK_RIGHT = 0x27,
		VK_DOWN = 0x28,
		VK_SELECT = 0x29,
		VK_PRINT = 0x2A,
		VK_EXECUTE = 0x2B,
		VK_SNAPSHOT = 0x2C,
		VK_INSERT = 0x2D,
		VK_DELETE = 0x2E,
		VK_HELP = 0x2F,
		//
		VK_LWIN = 0x5B,
		VK_RWIN = 0x5C,
		VK_APPS = 0x5D,
		//
		VK_SLEEP = 0x5F,
		//
		VK_NUMPAD0 = 0x60,
		VK_NUMPAD1 = 0x61,
		VK_NUMPAD2 = 0x62,
		VK_NUMPAD3 = 0x63,
		VK_NUMPAD4 = 0x64,
		VK_NUMPAD5 = 0x65,
		VK_NUMPAD6 = 0x66,
		VK_NUMPAD7 = 0x67,
		VK_NUMPAD8 = 0x68,
		VK_NUMPAD9 = 0x69,
		VK_MULTIPLY = 0x6A,
		VK_ADD = 0x6B,
		VK_SEPARATOR = 0x6C,
		VK_SUBTRACT = 0x6D,
		VK_DECIMAL = 0x6E,
		VK_DIVIDE = 0x6F,
		VK_F1 = 0x70,
		VK_F2 = 0x71,
		VK_F3 = 0x72,
		VK_F4 = 0x73,
		VK_F5 = 0x74,
		VK_F6 = 0x75,
		VK_F7 = 0x76,
		VK_F8 = 0x77,
		VK_F9 = 0x78,
		VK_F10 = 0x79,
		VK_F11 = 0x7A,
		VK_F12 = 0x7B,
		VK_F13 = 0x7C,
		VK_F14 = 0x7D,
		VK_F15 = 0x7E,
		VK_F16 = 0x7F,
		VK_F17 = 0x80,
		VK_F18 = 0x81,
		VK_F19 = 0x82,
		VK_F20 = 0x83,
		VK_F21 = 0x84,
		VK_F22 = 0x85,
		VK_F23 = 0x86,
		VK_F24 = 0x87,
		//
		VK_NUMLOCK = 0x90,
		VK_SCROLL = 0x91,
		//
		VK_OEM_NEC_EQUAL = 0x92,   // '=' key on numpad
								   //
		VK_OEM_FJ_JISHO = 0x92,   // 'Dictionary' key
		VK_OEM_FJ_MASSHOU = 0x93,   // 'Unregister word' key
		VK_OEM_FJ_TOUROKU = 0x94,   // 'Register word' key
		VK_OEM_FJ_LOYA = 0x95,   // 'Left OYAYUBI' key
		VK_OEM_FJ_ROYA = 0x96,   // 'Right OYAYUBI' key
								 //
		VK_LSHIFT = 0xA0,
		VK_RSHIFT = 0xA1,
		VK_LCONTROL = 0xA2,
		VK_RCONTROL = 0xA3,
		VK_LMENU = 0xA4,
		VK_RMENU = 0xA5,
		//
		VK_BROWSER_BACK = 0xA6,
		VK_BROWSER_FORWARD = 0xA7,
		VK_BROWSER_REFRESH = 0xA8,
		VK_BROWSER_STOP = 0xA9,
		VK_BROWSER_SEARCH = 0xAA,
		VK_BROWSER_FAVORITES = 0xAB,
		VK_BROWSER_HOME = 0xAC,
		//
		VK_VOLUME_MUTE = 0xAD,
		VK_VOLUME_DOWN = 0xAE,
		VK_VOLUME_UP = 0xAF,
		VK_MEDIA_NEXT_TRACK = 0xB0,
		VK_MEDIA_PREV_TRACK = 0xB1,
		VK_MEDIA_STOP = 0xB2,
		VK_MEDIA_PLAY_PAUSE = 0xB3,
		VK_LAUNCH_MAIL = 0xB4,
		VK_LAUNCH_MEDIA_SELECT = 0xB5,
		VK_LAUNCH_APP1 = 0xB6,
		VK_LAUNCH_APP2 = 0xB7,
		//
		VK_OEM_1 = 0xBA,   // ';:' for US
		VK_OEM_PLUS = 0xBB,   // '+' any country
		VK_OEM_COMMA = 0xBC,   // ',' any country
		VK_OEM_MINUS = 0xBD,   // '-' any country
		VK_OEM_PERIOD = 0xBE,   // '.' any country
		VK_OEM_2 = 0xBF,   // '/?' for US
		VK_OEM_3 = 0xC0,   // '`~' for US
						   //
		VK_OEM_4 = 0xDB,  //  '[{' for US
		VK_OEM_5 = 0xDC,  //  '\|' for US
		VK_OEM_6 = 0xDD,  //  ']}' for US
		VK_OEM_7 = 0xDE,  //  ''"' for US
		VK_OEM_8 = 0xDF,
		//
		VK_OEM_AX = 0xE1,  //  'AX' key on Japanese AX kbd
		VK_OEM_102 = 0xE2,  //  "<>" or "\|" on RT 102-key kbd.
		VK_ICO_HELP = 0xE3,  //  Help key on ICO
		VK_ICO_00 = 0xE4,  //  00 key on ICO
						   //
		VK_PROCESSKEY = 0xE5,
		//
		VK_ICO_CLEAR = 0xE6,
		//
		VK_PACKET = 0xE7,
		//
		VK_OEM_RESET = 0xE9,
		VK_OEM_JUMP = 0xEA,
		VK_OEM_PA1 = 0xEB,
		VK_OEM_PA2 = 0xEC,
		VK_OEM_PA3 = 0xED,
		VK_OEM_WSCTRL = 0xEE,
		VK_OEM_CUSEL = 0xEF,
		VK_OEM_ATTN = 0xF0,
		VK_OEM_FINISH = 0xF1,
		VK_OEM_COPY = 0xF2,
		VK_OEM_AUTO = 0xF3,
		VK_OEM_ENLW = 0xF4,
		VK_OEM_BACKTAB = 0xF5,
		//
		VK_ATTN = 0xF6,
		VK_CRSEL = 0xF7,
		VK_EXSEL = 0xF8,
		VK_EREOF = 0xF9,
		VK_PLAY = 0xFA,
		VK_ZOOM = 0xFB,
		VK_NONAME = 0xFC,
		VK_PA1 = 0xFD,
		VK_OEM_CLEAR = 0xFE
	}

	enum MouseEvent : int {
		MOUSEEVENTF_ABSOLUTE = 0x8000,  //Нормированные абсолютные координаты
		MOUSEEVENTF_LEFTDOWN = 0x0002,  //Нажатие на левую кнопку мыши
		MOUSEEVENTF_LEFTUP = 0x0004,    //Поднятие левой кнопки мыши
		MOUSEEVENTF_MOVE = 0x0001,      //перемещение указателя мыши
		MOUSEEVENTF_RIGHTDOWN = 0x0008, //Нажатие на правую кнопку мыши
		MOUSEEVENTF_RIGHTUP = 0x0010    //Поднятие правой кнопки мыши
	}

	enum MyMouseEvent : int {
		Mouse_Right_Up = 0,  //Нормированные абсолютные координаты
		Mouse_Right_Down = 1,  //Нажатие на левую кнопку мыши
		Mouse_Left_Up = 2,    //Поднятие левой кнопки мыши
		Mouse_Left_Down = 3,      //перемещение указателя мыши
		Mouse_Left_Click = 4, //Нажатие на правую кнопку мыши
		Mouse_Right_Click = 5    //Поднятие правой кнопки мыши
	}

	class CppLib {
		public static List<IntPtr> ChildHendel;
		const int BM_CLICK = 0x00F5;

		public const int WH_KEYBOARD_LL = 13;
		public delegate IntPtr LowLevelKeyboardProcDelegate(int nCode, IntPtr wParam, IntPtr lParam);
		public static IntPtr m_hHook;

		public delegate bool EnumWindowsProc(IntPtr hwnd, IntPtr lParam);
		delegate bool EnumWindowsCallback(IntPtr currentWindowHandle, ref SearchData searchData);

		[StructLayout(LayoutKind.Sequential)]
		public struct Rect {
			public int left;
			public int top;
			public int right;
			public int bottom;
		}

		struct SearchData {
			public string WindowText;
			public IntPtr ParentHandle;
			public IntPtr ResultHandle;
		}

		public struct ChildEll {
			public IntPtr Handle;
			public string Name;
			public string ClassEll;
		}

		public struct ProcessInfo {
			public IntPtr Handle;
			public string Name;
			public Rect Border;
		}

		[StructLayout(LayoutKind.Sequential)]
		public struct KeyboardHookStruct {
			public readonly int VirtualKeyCode;
			public readonly int ScanCode;
			public readonly int Flags;
			public readonly int Time;
			public readonly IntPtr ExtraInfo;
		}

		[DllImport("user32.dll", CharSet = CharSet.Unicode)]
		public static extern bool GetWindowRect(IntPtr hwnd, ref Rect rectangle);

		[DllImport("user32.dll", SetLastError = true)]
		public static extern IntPtr FindWindow(string lpClassName, string lpWindowName);

		[DllImport("user32.dll", SetLastError = true)]
		public static extern IntPtr FindWindowEx(IntPtr parentHandle, IntPtr childAfter, string className, string windowTitle);

		[DllImport("user32.dll", CharSet = CharSet.Auto)]
		public static extern int GetWindowTextLength(HandleRef hWnd);

		[DllImport("user32.dll", CharSet = CharSet.Auto)]
		public static extern int GetWindowText(HandleRef hWnd, StringBuilder lpString, int nMaxCount);

		[DllImport("user32.dll", SetLastError = true, CharSet = CharSet.Auto)]
		public static extern bool SetWindowText(IntPtr hwnd, String lpString);

		[DllImport("user32.dll")]
		static extern bool EnumWindows(EnumWindowsCallback callback, ref SearchData searchData);

		[DllImport("user32.dll")]
		public static extern IntPtr GetParent(IntPtr childHandle);

		[DllImport("user32.dll")]
		public static extern void GetWindowText(IntPtr handle, StringBuilder resultWindowText, int maxTextCapacity);

		[DllImport("user32.dll")]
		[return: MarshalAs(UnmanagedType.Bool)]
		public static extern bool EnumChildWindows(IntPtr hwndParent, EnumWindowsProc lpEnumFunc, IntPtr lParam);

		[DllImport("user32.dll", CharSet = CharSet.Auto)]
		static extern IntPtr SendMessage(IntPtr hWnd, int Msg, IntPtr wParam, [Out] IntPtr lParam);

		[DllImport("User32.dll")]
		public static extern IntPtr SendMessage(IntPtr hWnd, int Msg, int wParam, string lParam);

		[DllImport("user32.dll", CharSet = CharSet.Auto)]
		static extern IntPtr SendMessage(IntPtr hWnd, UInt32 Msg, IntPtr wParam, [Out] StringBuilder lParam);

		[DllImport("user32.dll")]
		public static extern int GetClassName(IntPtr hWnd, StringBuilder buf, int nMaxCount);

		[DllImport("user32.dll", SetLastError = true)]
		[return: MarshalAs(UnmanagedType.Bool)]
		public static extern bool IsWindowVisible(IntPtr hWnd);

		[System.Runtime.InteropServices.DllImport("user32.dll", CharSet = System.Runtime.InteropServices.CharSet.Auto,
		CallingConvention = System.Runtime.InteropServices.CallingConvention.StdCall)]
		public static extern void mouse_event(uint dwFlags, int dx, int dy, int dwData, int dwExtraInfo);

		[return: MarshalAs(UnmanagedType.Bool)]
		[DllImport("user32.dll", SetLastError = true, CharSet = CharSet.Auto)]
		static extern bool PostMessage(IntPtr hWnd, uint Msg, int wParam, int lParam);

		[DllImport("user32.dll", SetLastError = true)]
		public static extern IntPtr SetWindowsHookEx(int idHook, LowLevelKeyboardProcDelegate lpfn, IntPtr hMod, int dwThreadId);

		[DllImport("user32.dll", SetLastError = true)]
		public static extern bool UnhookWindowsHookEx(IntPtr hhk);

		[DllImport("Kernel32.dll", SetLastError = true)]
		public static extern IntPtr GetModuleHandle(IntPtr lpModuleName);

		[DllImport("user32.dll", SetLastError = true)]
		public static extern IntPtr CallNextHookEx(IntPtr hhk, int nCode, IntPtr wParam, IntPtr lParam);

		static bool Callback2(IntPtr hwnd, IntPtr lParam) {
			ChildHendel.Add(hwnd);
			return true;
		}

		static bool Callback(IntPtr currentWindowHandle, ref SearchData searchData) {
			bool continueEnumeration = true;

			IntPtr currentWindowParentHandle = GetParent(currentWindowHandle);

			if (currentWindowParentHandle == searchData.ParentHandle) {
				var windowText = new StringBuilder(1024);
				GetWindowText(currentWindowHandle, windowText, windowText.Capacity);

				if (windowText.ToString() == searchData.WindowText) {
					searchData.ResultHandle = currentWindowHandle;
					continueEnumeration = false;
				}
			}
			return continueEnumeration;
		}

		public static IntPtr GetChildWindowHandle(string windowText, IntPtr parentHandle) {
			var searchData = new SearchData { ParentHandle = parentHandle, WindowText = windowText };
			EnumWindows(Callback, ref searchData);
			return searchData.ResultHandle;
		}
		
		public static List<ChildEll> GetChildEllHandle(IntPtr hwndParent) {
			List<ChildEll> res = new List<ChildEll>();
			ChildHendel = new List<IntPtr>();
			EnumChildWindows(hwndParent, Callback2, IntPtr.Zero);
			for(int i = 0;i < ChildHendel.Count;i++) {
				int length = (int)SendMessage(ChildHendel[i], (int)WindowsMessages.WM_GETTEXTLENGTH, IntPtr.Zero, IntPtr.Zero);
				StringBuilder sb;
				if (length != 0) {
					sb = new StringBuilder(length + 1);
					SendMessage(ChildHendel[i], (int)WindowsMessages.WM_GETTEXT, (IntPtr)sb.Capacity, sb);
				}
				else {
					sb = new StringBuilder("");
				}
				StringBuilder strClassName = new StringBuilder(256);
				GetClassName(ChildHendel[i], strClassName, strClassName.Capacity);
				res.Add(new ChildEll { Handle = ChildHendel[i], Name = sb.ToString(), ClassEll = strClassName.ToString() });
			}
			return res;
		}

		public static List<ProcessInfo> GetProcessList() {
			List<ProcessInfo> res = new List<ProcessInfo>();
			Process[] processlist = Process.GetProcesses();
			foreach (Process process in processlist) {
				if (!String.IsNullOrEmpty(process.MainWindowTitle)) {
					Rect buff = new Rect();
					GetWindowRect(process.MainWindowHandle, ref buff);
					res.Add(new ProcessInfo { Handle = process.MainWindowHandle, Name = process.MainWindowTitle, Border = buff });
				}
				else if (process.ProcessName == "explorer") {
					Rect buff = new Rect();
					GetWindowRect(process.MainWindowHandle, ref buff);
					res.Add(new ProcessInfo { Handle = process.MainWindowHandle, Name = "Рабочий стол", Border = buff });
				}
			}
			return res;
		}

		public static void SetTitleText(IntPtr Handle, String Text) {
			SendMessage(Handle, (int)WindowsMessages.WM_SETTEXT, Text.Length, Text);
			PostMessage(Handle, (int)WindowsMessages.WM_KEYDOWN, (int)VirtualKeyStates.VK_END, 0);
			PostMessage(Handle, (int)WindowsMessages.WM_KEYDOWN, (int)VirtualKeyStates.VK_SPACE, 0);
		}

		public static String GetTitleText(IntPtr Handle) {
			int length = (int)SendMessage(Handle, (int)WindowsMessages.WM_GETTEXTLENGTH, IntPtr.Zero, IntPtr.Zero);
			StringBuilder sb = new StringBuilder(length + 1);
			SendMessage(Handle, (int)WindowsMessages.WM_GETTEXT, (IntPtr)sb.Capacity, sb);
			return sb.ToString();
		}

		public static void ClickButton(IntPtr Handle) {
			SendMessage(Handle, BM_CLICK, IntPtr.Zero, IntPtr.Zero);
		}

		public static List<ChildEll> GetNeedChildEll(IntPtr Perent, String[] Classes) {
			List<ChildEll> Childs = GetChildEllHandle(Perent);
			List<ChildEll> res = new List<ChildEll>();
			for(int i = 0;i <Childs.Count;i++) {
				for (int j = 0; j < Classes.Length; j++) {
					if (Classes[j] == Childs[i].ClassEll) {
						res.Add(Childs[i]);
						break;
					}
				}
			}
			return res;
		}

		public static void MouseClickOnDesktop(MyMouseEvent e) {
			uint MouseE = 0;
			switch (e) {
				case MyMouseEvent.Mouse_Right_Click:
					MouseE = (int)MouseEvent.MOUSEEVENTF_ABSOLUTE | (int)MouseEvent.MOUSEEVENTF_RIGHTDOWN | (int)MouseEvent.MOUSEEVENTF_RIGHTUP;
					break;
				case MyMouseEvent.Mouse_Left_Click:
					MouseE = (int)MouseEvent.MOUSEEVENTF_ABSOLUTE | (int)MouseEvent.MOUSEEVENTF_LEFTDOWN | (int)MouseEvent.MOUSEEVENTF_LEFTUP;
					break;
					case MyMouseEvent.Mouse_Right_Down:
					MouseE = (int)MouseEvent.MOUSEEVENTF_ABSOLUTE | (int)MouseEvent.MOUSEEVENTF_RIGHTDOWN;
					break;
				case MyMouseEvent.Mouse_Right_Up:
					MouseE = (int)MouseEvent.MOUSEEVENTF_ABSOLUTE | (int)MouseEvent.MOUSEEVENTF_RIGHTUP;
					break;
				case MyMouseEvent.Mouse_Left_Down:
					MouseE = (int)MouseEvent.MOUSEEVENTF_ABSOLUTE | (int)MouseEvent.MOUSEEVENTF_LEFTDOWN;
					break;
				case MyMouseEvent.Mouse_Left_Up:
					MouseE = (int)MouseEvent.MOUSEEVENTF_ABSOLUTE | (int)MouseEvent.MOUSEEVENTF_LEFTUP;
					break;
			}
			mouse_event(MouseE, System.Windows.Forms.Cursor.Position.X, System.Windows.Forms.Cursor.Position.Y, 0, 0);
		}

		public static void SetHook(LowLevelKeyboardProcDelegate m_callback) {
			m_hHook = SetWindowsHookEx(WH_KEYBOARD_LL,m_callback, GetModuleHandle(IntPtr.Zero), 0);
		}

		public static void Unhook() {
			UnhookWindowsHookEx(m_hHook);
		}
	}
}