package main

import (
	"fmt"
	"net"
	"os"
	"runtime"
)

func echo(conn net.Conn) {
	buf := make([]byte, 1024)
	for {
		n, err := conn.Read(buf)
                if err != nil {
                        break;
                }

                conn.Write(buf[:n])
                if err != nil {
                        fmt.Println("ERROR", err)
                        os.Exit(1)
                }
	}
}

func main() {
	runtime.GOMAXPROCS(1)

	l, err := net.Listen("tcp", "0.0.0.0:55555")
	if err != nil {
		fmt.Println("ERROR", err)
		os.Exit(1)
	}

	for {
		conn, err := l.Accept()
		if err != nil {
			fmt.Println("ERROR", err)
			continue
		}
		go echo(conn)
	}
}
